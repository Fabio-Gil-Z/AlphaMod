"""Functions to build a structure from a sequence in the alignment file.
   Usually called by setting AutoModel.generate_method."""

from modeller.alignment import Alignment
from modeller.selection import Selection


def _set_seq_id(mdl, aln):
    """Set the model's seq_id field using the most similar template.
       This is done automatically by transfer_xyz, but generate methods
       that don't use transfer_xyz to create the model from the templates
       will need to do this manually."""
    target = aln[-1]
    seq_id = 0.
    for known in aln[:-1]:
        seq_id = max(seq_id, known.get_sequence_identity(target))
    mdl.seq_id = seq_id


def generate_xyz(mdl, aln):
    """Build coordinates from scratch (ignore templates)"""

    mdl.read_top_par()
    mdl.create_topology(aln)
    mdl.build(initialize_xyz=True, build_method='3D_INTERPOLATION')
    mdl.make_valid_pdb_coordinates()
    _set_seq_id(mdl, aln)


def transfer_xyz(mdl, aln):
    """Build structure by copying equivalent coordinates from the templates"""

    # If initial malign3d was requested, orient the template structures but
    # then restore the original alignment
    if mdl.initial_malign3d:
        aln.clear()
        aln.append(file=mdl.alnfile, align_codes=mdl.knowns)
        aln.malign3d(fit=False, gap_penalties_3d=(0, 4))
        mdl.read_alignment(aln)
    mdl.read_top_par()
    mdl.create_topology(aln)
    mdl.transfer_xyz(aln, cluster_cut=-1.0)
    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    _rebuild_bad_rings(mdl)
    mdl.make_valid_pdb_coordinates()


def _detect_bad_ring(r, ring_atoms):
    """Return True if the given residue contains a bad phenyl ring.
       Given an ordered list of the 6 atoms in the ring, measure three atom
       pair distances across the ring. In an ideal structure, these should each
       be around 2.6-2.7 angstroms. If they're substantially less than this,
       mark the ring as 'bad'."""
    for i in range(3):
        try:
            a1 = r.atoms[ring_atoms[i]]
            a2 = r.atoms[ring_atoms[i + 3]]
        except KeyError:
            continue
        dx = a1.x - a2.x
        dy = a1.y - a2.y
        dz = a1.z - a2.z
        d = dx*dx + dy*dy + dz*dz
        if d < 5.0:
            return True


def _unbuild_ring(r, atoms):
    s = Selection()
    for a in atoms:
        try:
            s.add(r.atoms[a])
        except KeyError:
            continue
    s.unbuild()


def _rebuild_bad_rings(mdl):
    """Rebuild the structure of any distorted phenyl ring, from internal
       coordinates. This is done because the optimizer sometimes has a hard
       time recovering the ideal geometry if given badly distorted rings."""
    rings = {'PHE': (('CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2'),
                     ('HD1', 'HE1', 'HZ', 'HE2', 'HD2')),
             'TRP': (('CD2', 'CE3', 'CZ3', 'CH2', 'CZ2', 'CE2'),
                     ('HE3', 'HZ3', 'HH2', 'HZ2')),
             'TYR': (('CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2'),
                     ('HD1', 'HE1', 'OH', 'HH', 'HE2', 'HD2'))}
    to_rebuild = []
    for r in mdl.residues:
        if r.name in rings:
            ring_atoms, extra_atoms = rings[r.name]
            if _detect_bad_ring(r, ring_atoms):
                to_rebuild.append(r)
                _unbuild_ring(r, ring_atoms + extra_atoms)
    if len(to_rebuild) > 0:
        print("The following %d residues contain 6-membered rings with "
              "poor geometries\nafter transfer from templates. Rebuilding "
              "rings from internal coordinates:\n   %s"
              % (len(to_rebuild),
                 "\n   ".join([str(r) for r in to_rebuild])))
        mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')


def read_xyz(mdl, aln):
    """Read in the initial structure from an existing file"""

    # Create two copies of the template sequence:
    a = Alignment(mdl.env)
    a.append(file=mdl.alnfile, align_codes=[mdl.sequence]*2)

    # Use the initial model as the first structure in the alignment:
    a[0].prottyp = 'structureX'
    a[0].atom_file = a[0].code = mdl.my_inifile

    mdl.read_top_par()
    mdl.create_topology(aln)

    # Get model coordinates from the initial model, making sure that the
    # sequence is correct and any missing atoms are filled in:
    mdl.transfer_xyz(a, cluster_cut=-1.0)
    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

    # Note that transfer_xyz will set mdl.seq_id to 100%, which is not what
    # we want. So calculate it correctly here.
    _set_seq_id(mdl, aln)
