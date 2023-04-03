from modeller.util.modlist import AppendList
from modeller.symmetry import Symmetry
import _modeller


class SymmetryList(AppendList):

    def __init__(self, mdl):
        self.__mdl = mdl
        AppendList.__init__(self)

    def __len__(self):
        sym = _modeller.mod_model_sym_get(self.__mdl.modpt)
        return _modeller.mod_symmetry_nsegsym_get(sym)

    def _getfunc(self, indx):
        # temporary hack: return a null Symmetry object
        return Symmetry((), (), 1.0)

    def append(self, obj):
        if not isinstance(obj, Symmetry):
            raise TypeError("can only use Symmetry objects here")
        obj._add_segments(self.__mdl)

    def report(self, deviation):
        """Report symmetry restraint violations"""
        return _modeller.mod_symmetry_report(self.__mdl.modpt, deviation)
