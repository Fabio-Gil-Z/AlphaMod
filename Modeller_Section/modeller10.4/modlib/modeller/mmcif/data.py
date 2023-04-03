# Note that this code needs Python 2.6 or later, unlike much of Modeller
# which works as far back as 2.3, so imports need to be protected by a
# version check

import modeller.features


class _CifAtom(object):
    """Map a Modeller Atom to mmCIF ids"""
    def __init__(self, atom, entity_map):
        self._atom = atom
        chain = atom.residue.chain
        # Modeller chain index is 1-based
        self.entity_id = entity_map[chain.index - 1]
        # Get seq_id offset (index is from start of model, but we need it to
        # be from the start of the chain)
        self.chain_offset = chain.residues[0].index - 1
    # Evaluate only when needed
    atom_id = property(lambda self: self._atom.name)
    comp_id = property(lambda self: self._atom.residue.name)
    asym_id = property(lambda self: self._atom.residue.chain.name)
    seq_id = property(
        lambda self: self._atom.residue.index - self.chain_offset)


class TemplateChain(object):
    """A single chain of a template"""
    pass


class TargetChain(object):
    """A single chain in the target"""
    pass


class Alignment(object):
    """A multiple alignment between one TargetChain and one or
       more TemplateChains"""
    pass


class Data(object):
    """Some data used in the modeling."""
    pass


class AlignmentData(Data):
    """Alignment data used in the modeling."""
    content_type = 'alignment'


class RestraintData(Data):
    """Restraint data used in the modeling."""
    content_type = 'restraint'


class Restraints(object):
    """A set of restraints used in the modeling."""
    pass


class DistanceRestraints(Restraints, list):
    """A list of distance restraints used in the modeling."""
    restraint_type = 'Distance restraints'


class Protocol(list):
    """A list of ProtocolSteps"""
    pass


class ProtocolStep(object):
    """A step in a protocol"""
    software_used = 'MODELLER'


class ModelingProtocolStep(ProtocolStep):
    """A ProtocolStep where modeling is done"""
    method_type = 'modeling'


class CifData(object):
    def add_alignment_info(self, aln, knowns, seq, mdl):
        self.mdl = mdl
        self.data = []
        self._add_target_chains(aln[seq], mdl)
        self._add_template_chains(aln, knowns, aln[seq])
        self._add_alignments(aln)

    def add_restraints(self, user_restraints):
        self.restraints = []
        # Only handle single-feature restraints for now
        self._add_distance_restraints(
            r for r in user_restraints
            if isinstance(r._features, modeller.features.Distance))

    def add_modeling_protocol(self):
        # Right now we just have a single protocol, containing a single
        # step - for modeling
        p = Protocol()
        s = ModelingProtocolStep()
        s.id = 1
        s.input_data = self.data[:]
        p.id = 1
        p.append(s)
        self.protocols = [p]

    def _add_distance_restraints(self, rs):
        self.distance_restraints = DistanceRestraints(rs)
        if len(self.distance_restraints) > 0:
            self._add_restraint(self.distance_restraints)

    def _add_data(self, data):
        """Add some data used in the modeling."""
        self.data.append(data)
        data.id = len(self.data)

    def _add_restraint(self, rs):
        """Add a set of restraints."""
        self.restraints.append(rs)
        rs.id = len(self.restraints)

        # add Data for this set of restraints
        d = RestraintData()
        d.name = 'User-provided restraints'
        rs.data = d
        self._add_data(d)

    def _get_gapped_sequence(self, chain):
        """Return a sequence including gaps for the given chain"""
        # todo: this is not quite right because leading gaps for the N terminus
        # will include trailing gaps of the C terminus of the previous chain
        def rescode_with_leading_gap(residue):
            return '-' * residue.get_leading_gaps() + residue.code
        return ''.join(rescode_with_leading_gap(r) for r in chain.residues)

    def _add_alignments(self, aln):
        """Add an Alignment for each TargetChain."""
        # add Data for this alignment
        d = AlignmentData()
        d.name = 'Target Template Alignment'
        self._add_data(d)

        target_to_alignment = {}
        self.alignments = []
        for nchain, target in enumerate(self.target_chains):
            a = Alignment()
            target_to_alignment[target] = a
            a.data = d
            a.target = target
            a.templates = []
            a.id = nchain + 1
            self.alignments.append(a)
        # Add template information
        for template in self.template_chains:
            a = target_to_alignment[template.target_chain]
            a.templates.append(template)

    def _get_target_chain(self, chain, target):
        """Get the TargetChain object that aligns with this chain"""
        # We just return the first match. This will miss cases where a template
        # chain aligns with multiple target chains, but this isn't handled
        # by the dictionary anyway (and is usually a modeling error)
        for r in chain.residues:
            target_r = r.get_aligned_residue(target)
            if target_r:
                # Modeller chain index is 1-based
                return self.target_chains[target_r.chain.index - 1]

    def _add_template_chains(self, aln, knowns, target):
        ordinal = 1
        self.template_chains = []
        for k in knowns:
            seq = aln[k]
            for chain in seq.chains:
                t = TemplateChain()
                # Assume PDB code is first 4 characters of align code
                t.pdb_code = k[:4].upper()
                t.id = ordinal
                t.asym_id = chain.name
                # todo: handle insertion codes, 1-based numbering vs
                # PDB numbers?
                t.seq_range = (chain.residues[0].num,
                               chain.residues[-1].num)
                # todo: handle non-standard residues
                t.sequence = ''.join(r.code for r in chain.residues)
                t.gapped_sequence = self._get_gapped_sequence(chain)
                t.sequence_can = ''.join(r.code for r in chain.residues)
                t.target_chain = self._get_target_chain(chain, target)
                self.template_chains.append(t)
                ordinal += 1

    def _get_target_entities(self, seq):
        # Mapping from chain # to entity #
        self._entity_map = {}
        seen_seqs = {}
        for nchain, chain in enumerate(seq.chains):
            s = tuple(r.type for r in chain.residues)
            seen_seqs[s] = None
            self._entity_map[nchain] = len(seen_seqs)

    def _add_target_chains(self, seq, mdl):
        self._get_target_entities(seq)
        ordinal = 1
        self.target_chains = []
        for nchain, chain in enumerate(seq.chains):
            model_chain = mdl.chains[nchain]
            t = TargetChain()
            t.id = ordinal
            # Use output model chain IDs, not the original sequence (chain
            # IDs may have been changed)
            t.asym_id = model_chain.name
            t.entity_id = self._entity_map[nchain]
            t.sequence = ''.join(r.code for r in chain.residues)
            t.gapped_sequence = self._get_gapped_sequence(chain)
            # todo: check if this numbering is correct (should match that
            # in entity_poly_seq)
            t.seq_range = (1, len(chain.residues))
            self.target_chains.append(t)
            ordinal += 1

    def write_mmcif(self, writer):
        self._write_template_details(writer)
        self._write_target_details(writer)
        self._write_oligomer_modeling(writer)
        self._write_alignment(writer)
        self._write_data(writer)
        self._write_restraints(writer)
        self._write_distance_restraints(writer)
        self._write_protocols(writer)

    def _write_template_details(self, writer):
        with writer.loop('_ma_template_details',
                         ['id', 'group_id', 'target_entity_id',
                          'pdb_code', 'pdb_label_asym_id',
                          'pdb_seq_id_begin', 'pdb_seq_id_end',
                          'pdb_seq_one_letter_code',
                          'pdb_seq_one_letter_code_can']) as lp:
            for t in self.template_chains:
                # todo: group templates?
                lp.write(id=t.id, group_id=1,
                         target_entity_id=t.target_chain.entity_id,
                         pdb_code=t.pdb_code,
                         pdb_label_asym_id=t.asym_id,
                         pdb_seq_id_begin=t.seq_range[0],
                         pdb_seq_id_end=t.seq_range[1],
                         pdb_seq_one_letter_code=t.sequence,
                         pdb_seq_one_letter_code_can=t.sequence_can)

    def _write_target_details(self, writer):
        with writer.loop('_ma_target_details',
                         ['id', 'entity_id', 'seq_db_one_letter_code']) as lp:
            for t in self.target_chains:
                lp.write(id=t.id, entity_id=t.entity_id,
                         seq_db_one_letter_code=t.sequence)

    def _write_oligomer_modeling(self, writer):
        with writer.loop('_ma_target_oligomer_modeling',
                         ['ordinal_id', 'oligomer_group_id',
                          'target_entity_id', 'asym_id', 'seq_id_begin',
                          'seq_id_end']) as lp:
            ordinal = 1
            for t in self.target_chains:
                # Modeller always models a single model, so everything is in
                # the same oligomer group (id=1)
                lp.write(ordinal_id=ordinal, oligomer_group_id=1,
                         target_entity_id=t.entity_id, asym_id=t.asym_id,
                         seq_id_begin=t.seq_range[0],
                         seq_id_end=t.seq_range[1])
                ordinal += 1

    def _write_alignment(self, writer):
        # todo: populate with info on how the alignment was made
        with writer.loop('_ma_alignment_details',
                         ['id', 'alignment_id', 'target_id',
                          'template_id', 'data_id']) as lp:
            ordinal = 1
            for a in self.alignments:
                for template in a.templates:
                    lp.write(id=ordinal, alignment_id=a.id,
                             target_id=a.target.id,
                             template_id=template.id,
                             data_id=a.data.id)
                    ordinal += 1

        with writer.loop('_ma_alignment',
                         ['ordinal_id', 'alignment_id', 'target_template_flag',
                          'sequence']) as lp:
            ordinal = 1
            for a in self.alignments:
                for template in a.templates:
                    lp.write(ordinal_id=ordinal, alignment_id=a.id,
                             target_template_flag=2,  # Template
                             sequence=template.gapped_sequence)
                    ordinal += 1
                lp.write(ordinal_id=ordinal, alignment_id=a.id,
                         target_template_flag=1,  # Target
                         sequence=a.target.gapped_sequence)

    def _write_data(self, writer):
        with writer.loop('_ma_data',
                         ['id', 'name', 'content_type']) as lp:
            for d in self.data:
                lp.write(id=d.id, name=d.name, content_type=d.content_type)

    def _write_restraints(self, writer):
        with writer.loop('_ma_restraints',
                         ['ordinal_id', 'restraint_id', 'data_id', 'name',
                          'restraint_type', 'details']) as lp:
            ordinal = 1
            for r in self.restraints:
                lp.write(ordinal_id=ordinal, restraint_id=r.id,
                         data_id=r.data.id, restraint_type=r.restraint_type)
                ordinal += 1

    def _write_distance_restraints(self, writer):
        with writer.loop('_ma_distance_restraints',
                         ['ordinal_id', 'restraint_id', 'group_id',
                          'entity_id_1', 'asym_id_1', 'seq_id_1', 'comp_id_1',
                          'atom_id_1', 'entity_id_2', 'asym_id_2', 'seq_id_2',
                          'comp_id_2', 'atom_id_2',
                          'distance_threshold', 'uncertainty']) as lp:
            ordinal = 1
            for r in self.distance_restraints:
                # Assume all such restraints are Gaussian for now
                atoms = r._features.indices_to_atoms(
                    self.mdl, r._features.get_atom_indices()[0])
                atoms = [_CifAtom(a, self._entity_map) for a in atoms]
                lp.write(ordinal_id=ordinal,
                         restraint_id=self.distance_restraints.id,
                         group_id=1,
                         atom_id_1=atoms[0].atom_id,
                         comp_id_1=atoms[0].comp_id,
                         asym_id_1=atoms[0].asym_id,
                         seq_id_1=atoms[0].seq_id,
                         entity_id_1=atoms[0].entity_id,
                         atom_id_2=atoms[1].atom_id,
                         comp_id_2=atoms[1].comp_id,
                         asym_id_2=atoms[1].asym_id,
                         seq_id_2=atoms[1].seq_id,
                         entity_id_2=atoms[1].entity_id,
                         distance_threshold=r._parameters[0],
                         uncertainty=r._parameters[1])
                ordinal += 1

    def _write_protocols(self, writer):
        self._write_step_input(writer)
        self._write_protocol_steps(writer)

    def _write_step_input(self, writer):
        with writer.loop('_ma_step_input',
                         ['id', 'protocol_id', 'step_id', 'data_id',
                          'content_sub_type']) as lp:
            ordinal = 1
            for p in self.protocols:
                for s in p:
                    for d in s.input_data:
                        lp.write(id=ordinal, protocol_id=p.id, step_id=s.id,
                                 data_id=d.id)
                        ordinal += 1

    def _write_protocol_steps(self, writer):
        with writer.loop('_ma_protocol_step',
                         ['ordinal_id', 'protocol_id', 'step_id', 'asym_id',
                          'entity_id', 'method_type', 'step_name',
                          'software_used']) as lp:
            ordinal = 1
            for p in self.protocols:
                for s in p:
                    # Apply protocol step to all chains in the model:
                    for t in self.target_chains:
                        lp.write(ordinal_id=ordinal, protocol_id=p.id,
                                 step_id=s.id, asym_id=t.asym_id,
                                 entity_id=t.entity_id,
                                 method_type=s.method_type,
                                 software_used=s.software_used)
                        ordinal += 1
