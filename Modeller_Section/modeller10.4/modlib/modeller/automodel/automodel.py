"""Classes to build complete model(s) using template information"""

import modeller.restraints
from modeller.model import Model
from modeller.information import info
from modeller.alignment import Alignment
from modeller.error import ModellerError
from modeller.energy_data import EnergyData
from modeller.selection import Selection
import modeller.modfile as modfile
from modeller import log, physical, IOData
from modeller.scripts import align_strs_seq
from modeller.optimizers import ConjugateGradients, actions
from modeller.automodel import refine, generate, randomize, autosched
import sys
from modeller.util.deprecation import _deprecation_handler


_HAVE_IO = False
if sys.version_info[:2] >= (2, 6):
    try:
        import io
        import modeller.mmcif.data
        import modeller.mmcif.writer
        _HAVE_IO = True
    except ImportError:
        pass


class _RestraintsWrapper(modeller.restraints.Restraints):
    def __init__(self, mdl):
        super(_RestraintsWrapper, self).__init__(mdl)
        self._user_restraints = []

    def add(self, *args):
        self._user_restraints.extend(args)
        return super(_RestraintsWrapper, self).add(*args)


class AutoModel(Model):
    """Automatically build complete model(s) using template information"""

    restraints = None  # Override attribute in base class

    _alignment_info = None

    # Default variables
    max_ca_ca_distance = 14.0
    max_n_o_distance = 11.0
    max_sc_mc_distance = 5.5
    max_sc_sc_distance = 5.0
    create_restraints = True
    deviation = 4.0
    toplib = '${LIB}/top_heav.lib'
    parlib = '${LIB}/par.lib'
    spline_on_site = True
    initial_malign3d = False
    final_malign3d = False
    starting_model = 1
    ending_model = 1
    write_intermediates = False
    library_schedule = autosched.normal
    pdb_ext = '.pdb'
    output_model_format = 'PDB'
    repeat_optimization = 1
    fit_in_refine = 'NO_FIT'
    refine_hot_only = False
    rstrs_refined = 1
    max_molpdf = 1e7
    optimize_output = 'NO_REPORT'
    max_var_iterations = 200
    trace_output = 10
    accelrys = info.accelrys
    alnfile = ''
    knowns = None
    sequence = None
    root_name = None
    my_inifile = inifile = ''
    csrfile = ''
    schfile = ''
    generate_method = None
    rand_method = None
    md_level = None
    assess_methods = None
    outputs = None
    parallel_job = None
    library_restraints = None
    tracefile = None
    blank_single_chain = False

    def __init__(self, env, alnfile, knowns, sequence,
                 deviation=None, library_schedule=None, csrfile=None,
                 inifile=None, assess_methods=None, root_name=None):
        Model.__init__(self, env)
        self.alnfile = alnfile
        if isinstance(knowns, tuple):
            knowns = list(knowns)
        if not isinstance(knowns, list):
            knowns = [knowns]
        self.knowns = knowns
        self.sequence = sequence
        self.root_name = root_name or sequence
        self.assess_methods = assess_methods
        self.set_defaults()
        libs = self.env.libs
        if libs.topology.in_memory and libs.parameters.in_memory:
            log.warning('AutoModel',
  """Topology and/or parameter libraries already in memory. These will
                be used instead of the AutoModel defaults. If this is not what you
                want, clear them before creating the AutoModel object with
                env.libs.topology.clear() and env.libs.parameters.clear()""")
        if deviation:
            self.deviation = deviation
        if library_schedule:
            self.library_schedule = library_schedule
        if csrfile:
            self.csrfile = csrfile
            self.create_restraints = False
        if inifile:
            self.my_inifile = inifile
            self.generate_method = generate.read_xyz
        # Catch any user-added restraints at the Python level
        self.restraints = _RestraintsWrapper(self)

    def set_output_model_format(self, fmt):
        """Set the output format for models. Supported values are PDB
           and MMCIF (default PDB)."""
        if fmt == 'PDB':
            self.output_model_format = 'PDB'
            self.pdb_ext = '.pdb'
        elif fmt == 'MMCIF':
            self.output_model_format = 'MMCIF'
            self.pdb_ext = '.cif'
        else:
            raise ValueError("Unknown model format")

    def make(self, exit_stage=0):
        """Build all models"""

        self.outputs = []
        self.homcsr(exit_stage)
        # Exit early?
        if exit_stage >= 1:
            return
        # Read all restraints once for the whole job (except when loops are
        # done when restraints are read for each *.B9999???? model):
        self.rd_restraints()

        # getting model(s) (topology library must be in memory; ensured
        # now by one of the three GENERATE_METHOD routines):
        atmsel = self._check_select_atoms()
        self.multiple_models(atmsel)

        self.write_summary(self.outputs, 'models')

        if self.final_malign3d:
            self.fit_models_on_template()

    def set_defaults(self):
        """Set most default variables"""
        seqfile = modfile._sanitize_filename(self.root_name)
        self.inifile = seqfile + '.ini'
        self.csrfile = seqfile + '.rsr'
        self.schfile = seqfile + '.sch'
        self.generate_method = generate.transfer_xyz
        self.rand_method = randomize.xyz
        self.md_level = refine.very_fast

    def auto_align(self, matrix_file='family.mat', overhang=0,
                   write_fit=False):
        """Create an initial alignment for fully automated comparative
           modeling. Use only when you have high template-model sequence
           identity."""
        segfile = self.alnfile
        self.alnfile = segfile + '.ali'
        align_strs_seq(self.env, segfile, self.alnfile, self.knowns,
                       self.sequence, matrix_file, overhang, write_fit)

    def very_fast(self):
        """Call this routine before calling 'make()' if you want really fast
           optimization"""
        self.max_ca_ca_distance = 10.0
        self.max_n_o_distance = 6.0
        self.max_sc_mc_distance = 5.0
        self.max_sc_sc_distance = 4.5
        # Note that all models will be the same if you do not
        # change rand_method
        self.rand_method = None
        self.max_var_iterations = 50
        self.library_schedule = autosched.fastest
        self.md_level = None

    def write_summary(self, outputs, modeltyp):
        """Print out a summary of all generated models"""
        ok = []
        failed = []
        for mdl in outputs:
            if mdl['failure']:
                failed.append(mdl)
            else:
                ok.append(mdl)
        if ok:
            self.write_ok_summary(ok, modeltyp)
        if failed:
            self.write_failure_summary(failed, modeltyp)

    def write_ok_summary(self, all, modeltyp):
        """Print out a summary of all successfully generated models"""
        print("\n>> Summary of successfully produced %s:" % modeltyp)
        fields = [x for x in all[0].keys() if x.endswith(' score')]
        fields.sort()
        fields = ['molpdf'] + fields
        header = '%-25s ' % 'Filename' + " ".join(['%14s' % x for x in fields])
        print(header)
        print('-' * len(header))
        for mdl in all:
            text = '%-25s' % mdl['name']
            for field in fields:
                if isinstance(mdl[field], (tuple, list)):
                    text = text + ' %14.5f' % mdl[field][0]
                else:
                    text = text + ' %14.5f' % mdl[field]
            print(text)
        print('')

    def write_failure_summary(self, all, modeltyp):
        """Print out a summary of all failed models"""
        print("\n>> Summary of failed %s:" % modeltyp)
        for mdl in all:
            print("%-25s %s" % (mdl['name'], mdl['failure']))
        print('')

    def rd_restraints(self):
        """Read all restraints. You can override this in subclasses to read
           additional restraints."""
        self.restraints.clear()
        self.restraints.append(file=self.csrfile)

    def get_model_filename(self, root_name, id1, id2, file_ext):
        """Returns the model PDB name - usually of the form
           foo.B000X000Y.pdb"""
        return modfile.default(file_id='.B', file_ext=file_ext,
                               root_name=root_name, id1=id1, id2=id2)

    def use_parallel_job(self, job):
        """Split multiple model building across a parallel job"""
        self.parallel_job = job

    def use_library_restraints(self, librestraints):
        """Use a restraints library"""
        self.library_restraints = librestraints

    def multiple_models(self, atmsel):
        """Build all models, given all the previously generated restraints"""
        if self.parallel_job is not None:
            self.parallel_multiple_models(atmsel)
        else:
            for num in range(self.starting_model, self.ending_model + 1):
                self.outputs.append(self.single_model(atmsel, num))

    def __getstate__(self):
        d = Model.__getstate__(self)
        # Don't send job information over the network, as Python can't
        # pickle sockets
        if 'parallel_job' in d:
            del d['parallel_job']
        return d

    def parallel_multiple_models(self, atmsel):
        """Split the model building across all workers in a parallel job"""
        from modeller.automodel.parallel import ModelTask

        job = self.parallel_job
        for num in range(self.starting_model, self.ending_model + 1):
            job.queue_task(ModelTask(self, num, atmsel))
        self.outputs.extend(job.run_all_tasks())

    def read_initial_model(self):
        """Read the initial model from a file.
           Note that you are counting on some model arrays not being deleted by
           read() (ie the charge array generated by generate_topology())."""
        # Make sure we read every atom, since when we write out the model, it
        # writes out every atom, not just non-HET/non-hydrogen/non-water
        io = IOData(copy=self.env.io)
        io.hetatm = io.water = io.hydrogen = True
        self.read(file=self.inifile, io=io,
                  model_format=self.output_model_format, keep_disulfides=True)

    def randomize_initial_structure(self, atmsel):
        """Get and randomize the initial structure"""
        self.read_initial_model()
        if self.rand_method:
            self.rand_method(atmsel)

    def new_trace_file(self, num):
        """Open a new optimization trace file"""
        if self.trace_output > 0:
            filename = modfile.default(
                file_ext='', file_id='.D', root_name=self.root_name,
                id1=0, id2=num)
            return open(filename, 'w')
        else:
            return None

    def get_optimize_actions(self):
        """Get actions to carry out during the initial optimization.
           The default is to write to a trace file every trace_output steps."""
        act = []
        if self.trace_output > 0:
            act.append(actions.Trace(self.trace_output, self.tracefile))
        return act

    def get_refine_actions(self):
        """Get actions to carry out during refinement.
           The default is to write to a trace file every trace_output steps."""
        return self.get_optimize_actions()

    def single_model(self, atmsel, num, parallel=False):
        """Build a single optimized model from the initial model"""
        self.tracefile = self.new_trace_file(num)
        try:
            if parallel:
                self.read_top_par()
                self.create_topology(self.read_alignment())

            self.randomize_initial_structure(atmsel)

            if parallel:
                self.rd_restraints()

            if not hasattr(self.library_schedule, "make_for_model"):
                raise TypeError("library_schedule should now be a Schedule "
                                "object, not an integer as in older versions "
                                "of Modeller")
            sched = self.library_schedule.make_for_model(self)
            sched = sched * self.env.schedule_scale
            fh = open(self.schfile, "w")
            sched.write(fh)
            fh.close()

            self.write_int(0, num)
            filename = self.get_model_filename(self.root_name, 9999, num,
                                               self.pdb_ext)
            out = {'name': filename, 'num': num, 'failure': None}
            try:
                for irepeat in range(0, self.repeat_optimization):
                    self.single_model_pass(atmsel, num, sched)
                self.to_iupac()
            except (ModellerError, OverflowError):
                detail = sys.exc_info()[1]
                if len(str(detail)) > 0:
                    out['failure'] = detail
                else:
                    out['failure'] = 'Optimization failed'
                del detail
            else:
                self.model_analysis(atmsel, filename, out, num)
        finally:
            if self.tracefile:
                self.tracefile.close()
            del self.tracefile
        return out

    def model_analysis(self, atmsel, filename, out, num):
        """Energy evaluation and assessment, and write out the model"""
        if self.accelrys:
            # Write the final model (Accelrys wants it before calculating the
            # profiles, so that the Biso column contains the original
            # template-derived averages)
            self.write(file=filename, model_format=self.output_model_format)
            for (id, norm) in (('.E', False), ('.NE', True)):
                atmsel.energy(output='LONG ENERGY_PROFILE',
                              normalize_profile=norm,
                              file=modfile.default(file_id=id, file_ext='',
                                                   root_name=self.root_name,
                                                   id1=9999, id2=num))
            # The new request from Lisa/Azat to print out only
            # stereochemical restraint violations (6/24/03):
            # select only stereochemical restraints (maybe add dihedral
            # angles?):
            scal = physical.Values(default=0, bond=1, angle=1, dihedral=1,
                                   improper=1, soft_sphere=1,
                                   disulfide_distance=1, disulfide_angle=1,
                                   disulfide_dihedral=1)
            for (id, norm) in (('.ES', False), ('.NES', True)):
                e = atmsel.energy(output='ENERGY_PROFILE',
                                  normalize_profile=norm, schedule_scale=scal,
                                  file=modfile.default(
                                      file_id=id, file_ext='',
                                      root_name=self.root_name,
                                      id1=9999, id2=num))
            (out['molpdf'], out['pdfterms']) = e
            self.user_after_single_model()
            # Do model assessment if requested
            self.assess(atmsel, self.assess_methods, out)
        else:
            # Do model assessment if requested
            assess_keys = self.assess(atmsel, self.assess_methods, out)

            e = atmsel.energy(output='LONG VIOLATIONS_PROFILE',
                              file=modfile.default(file_id='.V', file_ext='',
                                                   root_name=self.root_name,
                                                   id1=9999, id2=num))

            (out['molpdf'], out['pdfterms']) = e
            self.user_after_single_model()

            # Write the final model; Biso contains the violations profile
            self.write_final_model(filename, assess_keys, out)

    def write_final_model(self, filename, assess_keys, out):
        """Write the final model to a file"""
        self.write(file=filename, model_format=self.output_model_format,
                   extra_data=self.get_extra_model_data(assess_keys, out))

    def _get_residue_from_index(self, ind):
        """Map an alignment index into a model residue"""
        r = self.residues[ind - 1]
        return "%s:%s" % (r.num, r.chain.name)

    if _HAVE_IO:
        def get_mmcif_model_data(self, assess_keys, out):
            """Get extra data to insert into mmCIF output files"""
            if sys.version_info[0] >= 3:
                sio = io.StringIO()
            else:
                sio = io.BytesIO()
            sio.write('\n')
            w = modeller.mmcif.writer.CifWriter(sio)
            d = modeller.mmcif.data.CifData()
            d.add_alignment_info(self.read_alignment(), self.knowns,
                                 self.sequence, self)
            d.add_restraints(self.restraints._user_restraints)
            d.add_modeling_protocol()
            d.write_mmcif(w)
            return sio.getvalue()
    else:
        def get_mmcif_model_data(self, assess_keys, out):
            return ''

    def get_extra_model_data(self, assess_keys, out):
        """Get extra data to insert into mmCIF or PDB output files"""
        # Use only templates that aligned with the model
        knowns = [x for x in self.knowns
                  if self._alignment_info[x] is not None]
        if self.output_model_format == 'MMCIF':
            def format_templates(knowns):
                fmts = []
                for known in knowns:
                    info, seq_id = self._alignment_info[known]
                    for i in range(0, len(info), 2):
                        fmts.append(
                            "%s %s %s %s %s %.3f"
                            % (known, info[i][0], info[i+1][0],
                               self._get_residue_from_index(info[i][1]),
                               self._get_residue_from_index(info[i+1][1]),
                               seq_id))
                return fmts
            templates = '\n#\nloop_\n_modeller_template.id\n' \
                        + '_modeller_template.name\n' \
                        + '_modeller_template.template_begin\n' \
                        + '_modeller_template.template_end\n' \
                        + '_modeller_template.target_begin\n' \
                        + '_modeller_template.target_end\n' \
                        + '_modeller_template.pct_seq_id\n' \
                        + "\n".join(["%d %s" % (x[0]+1, x[1]) for x
                                     in enumerate(format_templates(knowns))])
            if isinstance(self.alnfile, Alignment):
                alnfile = ''
            else:
                alnfile = '\n_modeller.alignment ' + self.alnfile
            script = '\n_modeller.script ' + sys.argv[0]
            scores = self.get_scores_model_data(assess_keys, out,
                                                "_modeller.%s %s", "_")
            return('_modeller.sequence ' + self.sequence + alnfile
                   + script + scores + templates
                   + self.get_mmcif_model_data(assess_keys, out))
        else:
            def format_templates(knowns):
                fmts = []
                for known in knowns:
                    info, seq_id = self._alignment_info[known]
                    for i in range(0, len(info), 2):
                        fmts.append(
                            "REMARK   6 TEMPLATE: %s %s - %s "
                            "MODELS %s - %s AT %.1f%%"
                            % (known, info[i][0], info[i+1][0],
                               self._get_residue_from_index(info[i][1]),
                               self._get_residue_from_index(info[i+1][1]),
                               seq_id))
                return fmts
            templates = "\n" + "\n".join(format_templates(knowns))
            if isinstance(self.alnfile, Alignment):
                alnfile = ''
            else:
                alnfile = '\nREMARK   6 ALIGNMENT: ' + self.alnfile
            script = '\nREMARK   6 SCRIPT: ' + sys.argv[0]
            scores = self.get_scores_model_data(assess_keys, out,
                                                "REMARK   6 %s: %s", " ")
            return('REMARK   6 SEQUENCE: ' + self.sequence + alnfile
                   + script + scores + templates)

    def get_scores_model_data(self, assess_keys, out, fmt, ws_replace,
                              spacer='\n'):
        """Get contribution to extra_data from assessment scores"""
        def get_score(score):
            if isinstance(score, (tuple, list)):
                score = score[0]
            return "%.5f" % score
        if assess_keys:
            return spacer + \
                   "\n".join([fmt % (k.replace(' ', ws_replace),
                                     get_score(out[k])) for k in assess_keys])
        else:
            return ''

    def assess(self, atmsel, methods, out=None):
        """Assess the model using all given methods"""
        assess_keys = []
        assess_list = methods
        if assess_list:
            if not isinstance(assess_list, (tuple, list)):
                assess_list = [assess_list]
            for method in assess_list:
                (key, value) = method(atmsel)
                assess_keys.append(key)
                if out is not None:
                    out[key] = value
        return assess_keys

    def single_model_pass(self, atmsel, num, sched):
        """Perform a single pass of model optimization"""
        actions = self.get_optimize_actions()
        for (numstep, step) in enumerate(sched):
            molpdf = step.optimize(atmsel, output=self.optimize_output,
                                   max_iterations=self.max_var_iterations,
                                   actions=actions)
            self.write_int(numstep + 1, num)
            # Also check for molpdf being NaN (depends on Python version;
            # on 2.3 x86 it evaluates as equal to everything; with 2.4 x86
            # it is not greater or smaller than anything)
            if molpdf > self.max_molpdf \
               or (molpdf == 0. and molpdf == 1.) \
               or (not molpdf >= 0. and not molpdf < 0):
                log.error('single_model',
                          "Obj. func. (%.3f) exceeded max_molpdf (%.3f) "
                          % (molpdf, self.max_molpdf))
        actions = self.get_refine_actions()
        self.refine(atmsel, actions)

    def write_int(self, id1, id2):
        """Write intermediate model file during optimization, if
           so requested"""
        if self.write_intermediates:
            self.write(file=self.get_model_filename(self.root_name, id1, id2,
                                                    self.pdb_ext),
                       model_format=self.output_model_format)

    def read_alignment(self, aln=None):
        """Read the template-sequence alignment needed for modeling"""
        codes = self.knowns+[self.sequence]
        if isinstance(self.alnfile, Alignment):
            aln_codes = [x.code for x in self.alnfile]
            if aln_codes != codes:
                raise ValueError("Alignment must have the sequences in "
                                 "this order ('knowns' followed by "
                                 "'sequence'): %s. They are actually %s"
                                 % (codes, aln_codes))
            return self.alnfile
        if aln is None:
            aln = Alignment(self.env)
        aln.clear()
        aln.append(file=self.alnfile, align_codes=codes)
        return aln

    def check_alignment(self, aln):
        """Check the alignment for sanity"""
        aln.check()

    def _get_sequence_info(self, aln, tgt, known):
        all_start_ends = []
        start_res = end_res = None

        for pos in aln.positions:
            target_res = pos.get_residue(tgt)
            known_res = pos.get_residue(known)
            if target_res and known_res:
                if end_res is not None and end_res[0].chain != known_res.chain:
                    all_start_ends.extend((start_res, end_res))
                    start_res = None
                end_res = (known_res, target_res)
                if start_res is None:
                    start_res = (known_res, target_res)
        if start_res and end_res:
            all_start_ends.extend((start_res, end_res))
        if all_start_ends:
            return ([("%s:%s" % (res[0].num, res[0].chain.name), res[1].index)
                     for res in all_start_ends],
                    known.get_sequence_identity(tgt))

    def _get_alignment_info(self, aln):
        """Get and store certain information about the alignment.
           This is added to the header of any output model."""
        self._alignment_info = {}
        for known in self.knowns:
            self._alignment_info[known] = self._get_sequence_info(
                aln, aln[self.sequence], aln[known])

    def homcsr(self, exit_stage):
        """Construct the initial model and restraints"""
        # Check the alignment
        aln = self.read_alignment()

        # Since in general we do not want to loose the original alignment file
        # (which is usually not a temporary scratch file):
        if self.accelrys:
            # Accelrys code here (Azat, you may want to add the .tmp.ali part
            # so that you do not change the input alignment file, unless you
            # want to have it changed here for some other use elsewhere):
            aln.write(file=self.alnfile)  # file='.tmp.ali'
            codes = [seq.code for seq in aln]
            aln.read(file=self.alnfile, align_codes=codes)
            # modfile.delete(file='.tmp.ali')

        self.check_alignment(aln)

        self._get_alignment_info(aln)

        # make topology and build/read the atom coordinates:
        self.make_initial_model(aln)

        # exit early?
        if exit_stage == 2:
            return

        # make and write the stereochemical, homology, and special restraints?
        if self.create_restraints:
            self.mkhomcsr(Selection(self), aln)
            self.restraints.condense()
            self.restraints.write(self.csrfile)

    def make_initial_model(self, aln):
        """Make initial model topology and build/read the atom coordinates"""
        self.generate_method(self, aln)
        self.write(file=self.inifile, model_format=self.output_model_format)
        self._check_model_hetatm_water()

    def _check_model_hetatm_water(self):
        """Check to see if env.io.hetatm or water are set if we're using
           HETATMs or waters"""
        if not self.env.io.hetatm or not self.env.io.water:
            water = Selection(self).only_water_residues()
            # Note that het is always a superset of water
            het = Selection(self).only_het_residues() - water
            warn = """You have at least one %s residue in your model, but
              IOData.%s is False. (This means that Modeller will not read
              any %s data from your templates, which is usually not what
              you want. To fix this, set env.io.%s = True before creating
              the AutoModel or LoopModel object.)"""
            if len(water) > 0 and not self.env.io.water:
                log.warning("_check_model_hetatm_water",
                            warn % ('water', 'water', 'water', 'water'))
            if len(het) > 0 and not self.env.io.hetatm:
                log.warning("_check_model_hetatm_water",
                            warn % ('HETATM', 'hetatm', 'HETATM', 'hetatm'))

    def build_charmm_restraints(self, atmsel, rsr, aln):
        """Build restraints from CHARMM libraries"""
        rsr.make(atmsel, restraint_type='stereo',
                 spline_on_site=self.spline_on_site,
                 residue_span_range=(0, 99999))

        rsr.make(atmsel, aln=aln, restraint_type='phi-psi_binormal',
                 spline_on_site=self.spline_on_site,
                 residue_span_range=(0, 99999))

        for type in ['omega', 'chi1', 'chi2', 'chi3', 'chi4']:
            rsr.make(atmsel, aln=aln, restraint_type=type+'_dihedral',
                     spline_range=4.0, spline_dx=0.3, spline_min_points=5,
                     spline_on_site=self.spline_on_site,
                     residue_span_range=(0, 99999))

    def build_library_restraints(self, atmsel, rsr, libraries):
        """Build restraints from restraint libraries"""
        if not isinstance(libraries, (tuple, list)):
            libraries = (libraries,)
        for lib in libraries:
            if hasattr(lib, '__call__'):
                func = lib
            else:
                func = lib.make_restraints
            func(atmsel, rsr, self.env.edat.nonbonded_sel_atoms)

    def _check_select_atoms(self):
        """Select atoms to be optimized, and check for sanity"""
        atmsel = self.select_atoms()
        if not hasattr(atmsel, "get_atom_indices"):
            raise ModellerError("you must return a Selection object " +
                                "from select_atoms")
        elif len(atmsel) == 0:
            raise ModellerError("no atoms selected for optimization")
        elif atmsel.get_model() is not self:
            raise ModellerError("selection is defined on the wrong model")
        elif len(atmsel) < len(self.atoms):
            print("%d (of %d total) atoms selected for optimization"
                  % (len(atmsel), len(self.atoms)))
        return atmsel

    def mkhomcsr(self, atmsel, aln):
        """Construct typical comparative modeling restraints"""
        rsr = self.restraints
        rsr.clear()
        if self.library_restraints is not None:
            self.build_library_restraints(atmsel, rsr, self.library_restraints)
        else:
            self.build_charmm_restraints(atmsel, rsr, aln)

        # Generate homology-derived distance restraints:
        self.distance_restraints(atmsel, aln)

        # Generate restraints on non standard residues:
        self.nonstd_restraints(aln)

        # Special restraints have to be called last so that possible
        # cis-proline changes are reflected in the current restraints:
        self.special_restraints(aln)

    def distance_restraints(self, atmsel, aln):
        """Construct homology-derived distance restraints"""
        rsr = self.restraints
        # Only do the standard residue types for CA, N, O, MNCH, SDCH dst rsrs
        # (no HET or BLK residue types):
        stdres = atmsel.only_std_residues()
        calpha = stdres.only_atom_types('CA')
        nitrogen = stdres.only_atom_types('N')
        oxygen = stdres.only_atom_types('O')
        mainchain = stdres.only_mainchain()
        sidechain = stdres - mainchain

        for (dmodel, maxdis, rsrrng, rsrsgn, rsrgrp, sel1, sel2, stdev) in \
            ((5, self.max_ca_ca_distance, (2, 99999), True,
              physical.ca_distance, calpha, calpha, (0, 1.0)),
             (6, self.max_n_o_distance,   (2, 99999), False,
              physical.n_o_distance, nitrogen, oxygen, (0, 1.0)),
             (6, self.max_sc_mc_distance, (1, 2), False,
              physical.sd_mn_distance, sidechain, mainchain, (0.5, 1.5)),
             (6, self.max_sc_sc_distance, (2, 99999), True,
              physical.sd_sd_distance, sidechain, sidechain, (0.5, 2.0))):
            if len(sel1) > 0 and len(sel2) > 0:
                rsr.make_distance(sel1, sel2, aln=aln,
                                  spline_on_site=self.spline_on_site,
                                  distance_rsr_model=dmodel,
                                  restraint_group=rsrgrp,
                                  maximal_distance=maxdis,
                                  residue_span_range=rsrrng,
                                  residue_span_sign=rsrsgn,
                                  restraint_stdev=stdev, spline_range=4.0,
                                  spline_dx=0.7, spline_min_points=5)

    def nonstd_restraints(self, aln):
        """Create restraints on HETATM and BLK residues."""
        # Select all HETATM residues plus any ATOM residues that have
        # no defined topology (generally speaking, BLK residues)
        allatoms = Selection(self)
        selhet = allatoms.only_het_residues() | allatoms.only_no_topology()

        rsrgrp = physical.xy_distance
        self.het_std_restraints(aln, selhet, 10.0, 2.3, rsrgrp)
        self.het_het_restraints(aln, selhet, 7.0, 2.3, rsrgrp)
        self.het_internal_restraints(aln, selhet, rsrgrp)

    def het_std_restraints(self, aln, selhet, ca_distance, bond_distance,
                           rsrgrp):
        """Constrain each atom in a HETATM or BLK residue by its distance to all
           protein atoms within ``bond_distance`` (these interactions are also
           removed from the nonbonded list) and also less strongly to all
           protein CA atoms that are within ``ca_distance``. The former
           maintains any protein-HETATM bonding, and the latter maintains
           the residue in roughly the right position."""

        rsr = self.restraints
        selstd = Selection(self).only_std_residues()
        selca = selstd.only_atom_types('CA')

        print("%d atoms in HETATM/BLK residues constrained\n" % len(selhet)
              + "to protein atoms within %.2f angstroms\n" % bond_distance
              + "and protein CA atoms within %.2f angstroms" % ca_distance)
        # Build the bonds first; this avoids duplicated CA-ligand bonds since
        # make_distance() will not build restraints that are already on the
        # nonbond exclusion list
        for (sel, dist, excl,
             strength) in ((selstd, bond_distance, bond_distance, 0.05),
                           (selca, ca_distance, 0.0, 0.2)):
            rsr.make_distance(sel, selhet, aln=aln, distance_rsr_model=7,
                              maximal_distance=dist,
                              spline_on_site=self.spline_on_site,
                              restraint_group=rsrgrp,
                              restraint_stdev=(strength, 0.0),
                              residue_span_range=(1, 99999),
                              residue_span_sign=False, exclude_distance=excl)

    def het_het_restraints(self, aln, selhet, maximal_distance, bond_distance,
                           rsrgrp):
        """Constrain atom-atom distances between different HETATM or BLK
           residues, to maintain their relative orientation. Distances less
           than ``bond_distance`` are assumed to be bonds, and so are
           restrained more strongly and also excluded from the
           nonbonded list."""

        rsr = self.restraints
        # Handle bonds first, to avoid duplicated restraints (the second call
        # to make_distance() will not build restraints that are already on the
        # nonbond exclusion list)
        for (dist, excl, strength) in ((bond_distance, bond_distance, 0.05),
                                       (maximal_distance, 0.0, 0.2)):
            rsr.make_distance(selhet, selhet, aln=aln, distance_rsr_model=7,
                              maximal_distance=dist,
                              spline_on_site=self.spline_on_site,
                              restraint_group=rsrgrp,
                              restraint_stdev=(strength, 0.0),
                              residue_span_range=(1, 99999),
                              residue_span_sign=True,
                              exclude_distance=excl)

    def het_internal_restraints(self, aln, selhet, rsrgrp):
        """Constrain internal distances within any HETATM or BLK
           residue that has no defined topology, to keep it rigid."""

        # Get all residues without defined topology:
        selhet = selhet.only_no_topology()

        # Intra-residue:
        print("%d atoms in residues without defined topology\n" % len(selhet)
              + "constrained to be rigid bodies")
        rsr = self.restraints
        rsr.make_distance(selhet, selhet, aln=aln,
                          distance_rsr_model=7, maximal_distance=10.0,
                          spline_on_site=self.spline_on_site,
                          restraint_group=rsrgrp, restraint_stdev=(0.05, 0.0),
                          residue_span_range=(0, 0), residue_span_sign=True)

    def special_restraints(self, aln):
        """This can be redefined by the user to add special restraints.
           In this class, it does nothing."""
        pass

    def special_patches(self, aln):
        """This can be redefined by the user to add additional patches
           (for example, user-defined disulfides). In this class, it
           does nothing."""
        pass

    def read_top_par(self):
        """Read in the topology and parameter libraries"""
        libs = self.env.libs
        if not (libs.topology.in_memory and libs.parameters.in_memory):
            libs.topology.read(file=self.toplib)
            libs.parameters.read(file=self.parlib)

    def create_topology(self, aln, sequence=None):
        """Build the topology for this model"""
        if sequence is None:
            sequence = self.sequence
        self.clear_topology()
        self.generate_topology(aln[sequence],
                               blank_single_chain=self.blank_single_chain)
        self.default_patches(aln)
        self.special_patches(aln)

    def default_patches(self, aln):
        """Derive the presence of disulfides from template structures (you can
           still define additional disulfides in the special_patches
           routine)"""
        self.patch_ss_templates(aln)

    def guess_atom_type(self, atom):
        """Assign a CHARMM type (used for soft sphere and Lennard Jones
           interactions) for the given atom. This is done here by guessing a
           suitable CHARMM type based on the atom name, but can be extended
           in subclasses. Usually called from guess_atom_types()."""
        full = {'MG': 'MG', 'ZN': 'ZN', 'FE': 'FE'}
        first = {'C': 'CT2', 'N':  'N', 'O': 'O', 'H': 'H', 'S': 'S', 'P': 'P'}
        chmtyp = full.get(atom.name, None) or first.get(atom.name[:1], None)
        if chmtyp:
            atom.type = chmtyp

    def guess_atom_types(self):
        """Assign a CHARMM type (used for soft sphere and Lennard Jones
           interactions) for all atoms in residues without defined topology
           (usually BLK residues) which would otherwise all use the UNDF
           atom type. Override guess_atom_type() to improve the guess used."""
        def fmt_typ(typ):
            if typ is None:
                return repr(typ)
            else:
                return typ.name
        sel = Selection(self).only_no_topology()
        old_types = []
        new_types = []
        for x in sel:
            old_types.append(x.type)
            self.guess_atom_type(x)
            new_types.append(x.type)
        if len(sel) > 0:
            print("\nThe following CHARMM atom type assignments were made:")
            print("      Atom                 Old type        New type")
            for atom, old_type, new_type in zip(sel, old_types, new_types):
                print("      %-20s %-15s %-15s"
                      % (str(atom), fmt_typ(old_type), fmt_typ(new_type)))

    def select_atoms(self):
        """Select atoms to be optimized in the model building procedure. By
           default, this selects all atoms, but you can redefine this routine
           to select a subset instead."""
        return Selection(self)

    def initial_refine_hot(self, atmsel):
        """Do some initial refinement of hotspots in the model"""
        viol_rc = physical.Values(default=999)
        stereo_typ = (physical.bond, physical.angle, physical.dihedral,
                      physical.improper, physical.disulfide_distance,
                      physical.disulfide_angle, physical.disulfide_dihedral)
        homol_typ = (physical.ca_distance, physical.n_o_distance,
                     physical.omega_dihedral, physical.sd_mn_distance,
                     physical.phi_psi_dihedral, physical.sd_sd_distance)
        if self.rstrs_refined == 0:
            # Refine only hotspots that have badly violated stereochemical
            # restraints:
            for typ in stereo_typ:
                viol_rc[typ] = 4
        elif self.rstrs_refined == 1:
            # Refine hotspots that have badly violated stereochemical
            # restraints and the important homology-derived restraints:
            for typ in stereo_typ + homol_typ:
                viol_rc[typ] = 4
        elif self.rstrs_refined == 2:
            # Refine hotspots that have badly violated any kind of
            # restraints
            viol_rc['default'] = 4

        # Pick hot atoms (must pick whole residues because of sidechains):
        atmsel = atmsel.hot_atoms(pick_hot_cutoff=4.5, viol_report_cut=viol_rc)
        atmsel = atmsel.by_residue()
        # Pick all corresponding (violated and others) restraints:
        self.restraints.unpick_all()
        self.restraints.pick(atmsel)

        # Local optimization to prevent MD explosions:
        cg = ConjugateGradients()
        cg.optimize(atmsel, max_iterations=100, output=self.optimize_output)
        return atmsel

    def final_refine_hot(self, hot_atmsel, atmsel):
        """Do some final refinement of hotspots in the model"""
        # Get conjugate gradients refined hot spots:
        cg = ConjugateGradients()
        cg.optimize(hot_atmsel, max_iterations=200,
                    output=self.optimize_output)

        # Get all static restraints again and select all atoms
        self.restraints.unpick_all()
        self.restraints.pick(atmsel)

    def refine(self, atmsel, actions):
        """Refine the optimized model with MD and CG"""
        # Save the current model:
        if self.fit_in_refine != 'NO_FIT':
            self.write(file='TO_BE_REFINED.TMP',
                       model_format=self.output_model_format)

        # Possibly skip selecting hot atoms only and optimize all atoms:
        if self.refine_hot_only:
            hot_atmsel = self.initial_refine_hot(atmsel)
        else:
            hot_atmsel = atmsel

        # Do simulated annealing MD:
        if self.md_level:
            self.md_level(hot_atmsel, actions)

        # Possibly skip 'HOT CG' after MD:
        if self.refine_hot_only:
            self.final_refine_hot(hot_atmsel, atmsel)

        self.finish_refine(atmsel, actions)

        if 'NO_FIT' not in self.fit_in_refine:
            self.fit_refined('TO_BE_REFINED.TMP')

    def fit_refined(self, unrefined):
        """Evaluate gross changes between the initial and final
           refined model"""
        aln = Alignment(self.env)
        mdl2 = Model(self.env, file=unrefined)
        sel = Selection(self).only_atom_types('CA')
        sel.superpose(mdl2, aln)
        sel = Selection(self)
        sel.superpose(mdl2, aln)
        modfile.delete(unrefined)

    def finish_refine(self, atmsel, actions):
        """Complete refinement by locally optimizing with conjugate
           gradients"""
        cg = ConjugateGradients()
        cg.optimize(atmsel, max_iterations=200, output=self.optimize_output,
                    actions=actions)

    def user_after_single_model(self):
        """Used for any user analysis after building each model. Redefine as you
           see fit."""
        pass

    def align_models(self, aln):
        """Adds the sequences of all generated models to the alignment."""
        for mdl in [a for a in self.outputs if a['failure'] is None]:
            code = "%s_9999%04d" % (self.sequence, mdl['num'])
            aln.append_model(mdl=self, align_codes=code,
                             atom_files=mdl['name'])

    def cluster(self, cluster_cut=1.5):
        """Cluster all output models, and output an optimized cluster
           average"""
        self.read_initial_model()
        aln = Alignment(self.env)
        self.align_models(aln)
        if len(aln) == 0:
            log.error('cluster', 'No generated models - nothing to cluster!')
        aln.malign3d(gap_penalties_3d=(0, 3), fit=False)
        aln.append_model(mdl=self, align_codes='cluster',
                         atom_files='cluster.opt')
        self.transfer_xyz(aln, cluster_cut=cluster_cut)
        self.write(file='cluster.ini', model_format=self.output_model_format)
        self.read_top_par()
        self.rd_restraints()
        self.create_topology(aln, sequence='cluster')
        atmsel = self._check_select_atoms()
        self.restraints.unpick_all()
        self.restraints.pick(atmsel)
        self.restraints.condense()
        edat = EnergyData(copy=self.env.edat)
        edat.nonbonded_sel_atoms = 1
        atmsel.energy(output='LONG', edat=edat)
        cg = ConjugateGradients()
        fh = open('cluster.deb', 'w')
        cg.optimize(atmsel, actions=actions.Trace(5, fh),
                    max_iterations=self.max_var_iterations)
        fh.close()
        atmsel.energy()
        self.write(file='cluster.opt', model_format=self.output_model_format)
        aln.compare_structures(fit=True)

    def fit_models_on_template(self):
        """Superpose each of the generated models on the templates"""
        aln = Alignment(self.env)
        aln.append(file=self.alnfile, align_codes=self.knowns)
        self.align_models(aln)
        # To take care of the '.' in segment specs:
        aln.write(file='.tmp.ali', alignment_format='PIR')
        codes = [seq.code for seq in aln]
        aln.read(file='.tmp.ali', alignment_format='PIR', align_codes=codes)
        modfile.delete('.tmp.ali')

        aln.compare_structures(fit=True, output='SHORT', fit_atoms='CA')
        aln.malign3d(gap_penalties_3d=(0, 3), write_whole_pdb=False,
                     write_fit=True, fit=False, fit_atoms='CA',
                     current_directory=True,
                     edit_file_ext=(self.pdb_ext, '_fit' + self.pdb_ext),
                     output='SHORT')


# Modeller 9 compatibility
class automodel(AutoModel):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(automodel)
        AutoModel.__init__(self, *args, **keys)
