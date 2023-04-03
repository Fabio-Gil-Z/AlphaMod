from modeller.util.modlist import VarList
from modeller.excluded_pair import ExcludedPair
import _modeller


class ExcludedPairList(VarList):
    def __init__(self, mdl):
        self.__mdl = mdl
        VarList.__init__(self)

    def __len__(self):
        return _modeller.mod_model_nexcl_get(self.__mdl.modpt)

    def _setdimfunc(self, num):
        _modeller.mod_model_nexcl_set(self.__mdl.modpt, num)

    def _getfunc(self, indx):
        mdl = self.__mdl
        inds = _modeller.mod_excluded_pair_get(mdl.modpt, indx)
        p = ExcludedPair(*[mdl.atoms[i-1] for i in inds])
        return p

    def _setfunc(self, indx, obj):
        mdl = self.__mdl
        if not isinstance(obj, ExcludedPair):
            raise TypeError("can only use ExcludedPair objects here")
        _modeller.mod_excluded_pair_set(mdl.modpt, indx,
                                        obj._get_base_atoms(mdl))
