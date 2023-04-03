"""Classes for handling databases of protein sequences"""

import _modeller
from modeller.util import modutil
from modeller import residue
from modeller.util.modobject import ModObject
from modeller.util.deprecation import _deprecation_handler


class SequenceDB(ModObject):
    """Holds a database of protein sequences"""
    __modpt = None
    __wndpt = None
    __free_func = _modeller.mod_sequence_db_free
    __window_free_func = _modeller.mod_sequence_db_window_free
    __close_func = _modeller.mod_sequence_db_close
    env = None
    window_size = 1

    def __init__(self, env, **vars):
        self.__modpt = _modeller.mod_sequence_db_new(self)
        self.env = env.copy()
        if len(vars) > 0:
            self.read(**vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_sequence_db_new(self)

    def __del__(self):
        if self.__modpt:
            self.close()
            self.__free_func(self.__modpt)

    def __close_window(self):
        if self.__wndpt:
            self.__window_free_func(self.__wndpt, self.__modpt)
            self.__wndpt = None

    def __repr__(self):
        if len(self) == 1:
            return "Database of 1 sequence"
        else:
            return "Database of %d sequences" % len(self)

    def __str__(self):
        return "<%s>" % repr(self)

    def __get_modpt(self):
        return self.__modpt

    def __get_wndpt(self):
        if self.__wndpt is None:
            self.__wndpt = _modeller.mod_sequence_db_window_new(self.__modpt)
        return self.__wndpt

    def __len__(self):
        return _modeller.mod_sequence_db_nchn_get(self.modpt)

    def close(self):
        """Close any currently-open database"""
        self.__close_window()
        self.__close_func(self.modpt)

    def read(self, chains_list, seq_database_file, seq_database_format,
             clean_sequences=True, minmax_db_seq_len=(0, 999999)):
        """Reads in a database from a file"""
        self.__close_window()
        f = _modeller.mod_sequence_db_read
        return f(self.modpt, self.env.libs.modpt, chains_list,
                 seq_database_file, seq_database_format, clean_sequences,
                 minmax_db_seq_len)

    def convert(self, chains_list, seq_database_file, seq_database_format,
                outfile, clean_sequences=True, minmax_db_seq_len=(0, 999999)):
        """Converts a database to binary format"""
        self.__close_window()
        f = _modeller.mod_sequence_db_convert
        return f(self.modpt, self.env.libs.modpt, chains_list,
                 seq_database_file, seq_database_format, outfile,
                 clean_sequences, minmax_db_seq_len)

    def write(self, chains_list, seq_database_file, seq_database_format,
              window_size=1024):
        """Writes out a database to a file"""
        return _modeller.mod_sequence_db_write(
            self.modpt, self.env.libs.modpt, chains_list, seq_database_file,
            seq_database_format, window_size)

    def search(self, aln, seq_database_file, search_group_list,
               search_randomizations=0, search_top_list=20, off_diagonal=100,
               overhang=0, gap_penalties_1d=(-900., -50.),
               signif_cutoff=(4.0, 5.0), rr_file='$(LIB)/as1.sim.mat',
               matrix_offset=0., fast_search_cutoff=1.0, data_file=False,
               search_sort='LONGER', output='LONG',
               alignment_features='INDICES CONSERVATION',
               local_alignment=False, fast_search=False, window_size=1024,
               io=None, **vars):
        """Search for similar sequences"""
        if io is None:
            io = self.env.io
        func = _modeller.mod_sequence_db_search
        return func(self.modpt, aln.modpt, io.modpt, self.env.libs.modpt,
                    search_randomizations, search_top_list, off_diagonal,
                    overhang, gap_penalties_1d, signif_cutoff, rr_file,
                    matrix_offset, fast_search_cutoff, data_file,
                    search_group_list, search_sort, output, alignment_features,
                    seq_database_file, local_alignment, fast_search,
                    window_size)

    def filter(self, seqid_cut, output_grp_file, output_cod_file,
               gap_penalties_1d=(-900., -50.), matrix_offset=0.,
               rr_file='$(LIB)/as1.sim.mat', max_diff_res=30, window_size=512):
        """Cluster sequences by sequence-identity"""
        return _modeller.mod_sequence_db_filter(
            self.modpt, self.env.libs.modpt, gap_penalties_1d, matrix_offset,
            rr_file, seqid_cut, max_diff_res, output_grp_file,
            output_cod_file, window_size)

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx)
        if isinstance(ret, int):
            return Sequence(self, indx)
        else:
            return [self[ind] for ind in ret]

    modpt = property(__get_modpt)
    _wndpt = property(__get_wndpt)


class Sequence(object):
    """A single sequence in the database"""

    def __init__(self, sdb, num):
        self.sdb = sdb
        self.env = sdb.env
        self.num = num

    def __len__(self):
        nseq = _modeller.mod_sequence_db_nseqdb_get(self.sdb.modpt)
        return _modeller.mod_int1_get(nseq, self.num)

    def __repr__(self):
        if len(self) == 1:
            return "Sequence of 1 residue"
        else:
            return "Sequence of %d residues" % len(self)

    def __str__(self):
        return "<%s>" % repr(self)

    def _get_chain(self, num):
        ichn = _modeller.mod_sequence_db_chain_get(num, self.sdb.modpt,
                                                   self.sdb._wndpt,
                                                   self.sdb.window_size)
        return self.sdb._wndpt, ichn

    def __get_code(self):
        wndpt, ichn = self._get_chain(self.num)
        return _modeller.mod_sequence_db_code_get(wndpt, ichn)

    def __get_prottyp(self):
        wndpt, ichn = self._get_chain(self.num)
        return _modeller.mod_sequence_db_prottyp_get(wndpt, ichn)

    def __get_resol(self):
        wndpt, ichn = self._get_chain(self.num)
        resoldb = _modeller.mod_sequence_db_window_resol_get(wndpt)
        return _modeller.mod_float1_get(resoldb, ichn)

    def __get_residues(self):
        return ResidueList(self)

    code = property(__get_code, doc="Alignment code")
    prottyp = property(__get_prottyp, doc="Protein type")
    resol = property(__get_resol, doc="Resolution")
    residues = property(__get_residues, doc="All residues in this sequence""")


class ResidueList(object):
    """A list of residues (:class:`Residue` objects) in a single sequence in
       a sequence database."""

    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx)
        if isinstance(ret, int):
            return Residue(self.seq, ret)
        else:
            return [self[ind] for ind in ret]


class Residue(residue.Residue):
    """A single residue in a sequence database."""

    def __get_type(self):
        wndpt, ichn = self.mdl._get_chain(self.mdl.num)
        return _modeller.mod_sequence_db_restype_get(wndpt, ichn, self._num)
    type = property(__get_type, doc="Integer residue type")


# Modeller 9 compatibility
class sequence_db(SequenceDB):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(sequence_db)
        SequenceDB.__init__(self, *args, **keys)
