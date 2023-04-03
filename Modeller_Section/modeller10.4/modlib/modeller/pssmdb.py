import _modeller
from modeller.util.modobject import ModObject
from modeller.util.deprecation import _deprecation_handler


class PSSMDB(ModObject):
    """Holds a database of PSSMs"""
    __modpt = None
    __free_func = _modeller.mod_pssmdbobj_free
    env = None

    def __init__(self, env, **vars):
        self.__modpt = _modeller.mod_pssmdbobj_new(self)
        self.env = env.copy()
        if len(vars) > 0:
            self.read(**vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_pssmdbobj_new(self)

    def __del__(self):
        if self.__modpt:
            self.__free_func(self.__modpt)

    def __get_modpt(self):
        return self.__modpt

    def read(self, pssmdb_name, pssmdb_format):
        return _modeller.mod_pssmdbobj_read(self.__modpt, pssmdb_name,
                                            pssmdb_format)

    modpt = property(__get_modpt)


# Modeller 9 compatibility
class pssmdb(PSSMDB):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(pssmdb)
        PSSMDB.__init__(self, *args, **keys)
