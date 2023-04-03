import modeller
import sys


class SalignData(object):
    """Data returned from the 'Alignment.salign' method"""

    def __init__(self, aln_score, qscorepct):
        self.aln_score = aln_score
        self.qscorepct = qscorepct


def _get_align_type(align_block):
    if align_block == 0:
        return 'tree'
    else:
        return 'pairwise'


def _salign_fw_local_gaps1(aln, feature_weights, ogp, egp, align_block,
                           matrix_offset):
    """Local alignment with given parameters"""
    return aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                      rr_file='$(LIB)/as1.sim.mat', overhang=0,
                      gap_penalties_1d=(ogp, egp),
                      align_block=align_block,
                      local_alignment=True, matrix_offset=matrix_offset,
                      matrix_offset_3d=-0.5, gap_penalties_3d=(0, 3),
                      gap_gap_score=0, gap_residue_score=0,
                      alignment_type=_get_align_type(align_block), nsegm=2,
                      feature_weights=feature_weights,
                      improve_alignment=True, fit=True, write_fit=False,
                      output='ALIGNMENT QUALITY')


def _salign_fw_gaps3(aln, feature_weights, ogp3d, egp3d, align_block):
    """Global alignment with given parameters"""
    ogp = ogp3d
    egp = egp3d

    def _salign_wrap(aln, **keys):
        try:
            return aln.salign(auto_overhang=True, overhang_auto_limit=5,
                              overhang_factor=1, **keys)
        except modeller.ModellerError:
            print("SALIGN with auto_overhang failed: %s"
                  % str(sys.exc_info()[1]))
            print("Retrying without auto_overhang")
            return aln.salign(**keys)

    return _salign_wrap(aln, rms_cutoff=3.5, normalize_pp_scores=False,
                        rr_file='$(LIB)/as1.sim.mat', overhang=0,
                        gap_penalties_1d=(ogp, egp),
                        align_block=align_block,
                        local_alignment=False, matrix_offset=-0.2,
                        gap_penalties_3d=(ogp3d, egp3d), gap_gap_score=0,
                        gap_residue_score=0,
                        alignment_type=_get_align_type(align_block),
                        nsegm=2, feature_weights=feature_weights,
                        improve_alignment=True, fit=True, write_fit=False,
                        write_whole_pdb=False, output='ALIGNMENT QUALITY')


def _frange(start, end=None, inc=1.0):
    """A range function that accepts floating point increments"""
    if end is None:
        end = float(start)
        start = 0.0
    else:
        start = float(start)
    count = int((end - start)/inc)
    if start + (count*inc) != end:
        count += 1
    for i in range(count):
        yield start + i*inc


class _TemporaryDirectory(object):
    """Create a temporary directory, and delete it when this object
       goes out of scope."""
    def __init__(self):
        import tempfile
        import shutil
        self.shutil = shutil
        self.tmpdir = tempfile.mkdtemp()

    def __del__(self):
        if hasattr(self, 'shutil'):
            self.shutil.rmtree(self.tmpdir)

    def get_path(self, path):
        """Return the name of a file in the temporary directory"""
        import os
        return os.path.join(self.tmpdir, path)


class _BestAlignment(object):
    def __init__(self, aln):
        self.qscore = 0.
        self.aln = aln
        self.found_struc_align = False

    def try_seq_align(self, input_aln_file, output_aln_file,
                      weights, open_penalty, extend_penalty, align_block,
                      matrix_offset):
        self._try_candidate(input_aln_file, output_aln_file,
                            _salign_fw_local_gaps1, weights,
                            open_penalty, extend_penalty, align_block,
                            matrix_offset)

    def try_struc_align(self, input_aln_file, output_aln_file,
                        weights, open_penalty, extend_penalty, align_block):
        if self._try_candidate(input_aln_file, output_aln_file,
                               _salign_fw_gaps3, weights,
                               open_penalty, extend_penalty, align_block):
            self.found_struc_align = True

    def _try_candidate(self, input_aln_file, output_aln_file, aln_func,
                       weights, open_penalty, extend_penalty, align_block,
                       *args):
        better = False
        self.aln.clear()
        self.aln.append(file=input_aln_file)
        try:
            res = aln_func(self.aln, weights, open_penalty, extend_penalty,
                           align_block, *args)
            if res.qscorepct >= self.qscore:
                self.qscore = res.qscorepct
                self.aln.write(file=output_aln_file, alignment_format='PIR')
                better = True
            print("Qlty scrs %g\t%g\t%g" % (open_penalty, extend_penalty,
                                            res.qscorepct))
        except modeller.ModellerError:
            print("Set of parameters %s %g %g resulted in the "
                  "following error\t%s" % (str(weights), open_penalty,
                                           extend_penalty,
                                           str(sys.exc_info()[1])))
        return better


def iterative_structural_align(aln, align_block=0):
    """Given an alignment of structures, iterate over parameter values
       to obtain the best structural alignment."""
    for seq in aln:
        if seq.code == '_fix_pos':
            raise modeller.ModellerError("This method does not work with fixed"
                                         " alignment positions (_fix_pos)")
        if not hasattr(seq, 'atoms'):
            raise modeller.ModellerError("This method only works for an "
                                         "alignment of structures.")
    tmpdir = _TemporaryDirectory()
    fil = tmpdir.get_path("inp.pir")
    aln.write(file=fil)

    opfile = tmpdir.get_path("salign_local_mid.ali")
    opfile2 = tmpdir.get_path("salign_local.ali")

    # -- Iterating over values of gap penalties and matrix offset
    fw1 = (1., 0., 0., 0., 1., 0.)
    fw2 = (0., 1., 0., 0., 0., 0.)
    fw3 = (0., 0., 0., 0., 1., 0.)
    best = _BestAlignment(aln)

    # -- Iterating over gap penalties 1D to get initial alignments
    print("Iterate over 1D penalties to get initial alignments")
    for ogp in _frange(-150, 1, 50):
        for egp in _frange(-50, 1, 50):
            for mo in _frange(-3.0, -0.05, 0.3):
                best.try_seq_align(fil, opfile, fw1, ogp, egp, align_block, mo)

    # -- Iterating over gap penalties 3D to get final alignments
    print("Iterate over 3D penalties to get final alignments")
    for ogp3d in _frange(0, 3, 1):
        for egp3d in range(2, 5, 1):
            best.try_struc_align(opfile, opfile2, fw2, ogp3d, egp3d,
                                 align_block)

    # try alternate initial alignments only if the qmax score is less than 70%
    if best.qscore <= 70:
        print("Trying alternate initial alignments")
        for ogp in _frange(0.0, 2.2, 0.3):
            for egp in _frange(0.1, 2.3, 0.3):
                for mo in _frange(-3.0, -0.05, 0.3):
                    best.try_seq_align(fil, opfile, fw3, ogp, egp, align_block,
                                       mo)

        # -- Iterating over gap penalties 3D to get final alignments
        print("Trying alternate final alignments")
        for ogp3d in _frange(0, 3, 1):
            for egp3d in range(2, 5, 1):
                best.try_struc_align(opfile, opfile2, fw2, ogp3d, egp3d,
                                     align_block)

    print("final max quality = %g" % best.qscore)

    if best.found_struc_align:
        aln.clear()
        aln.append(file=opfile2)
    else:
        raise modeller.ModellerError("Structure alignment failed")
