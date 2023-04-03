from modeller.util.deprecation import _deprecation_handler


class ModObject(object):
    def __setattr__(self, name, val):
        if name not in dir(self):
            print("runcmd_____W>: creation of new member '%s' in %s: "
                  "possible typo!" % (name, str(self)))
        object.__setattr__(self, name, val)

    def __getstate__(self):
        # Python's pickle won't work with internal Modeller objects,
        # so remove them from the list of objects to be pickled (don't
        # delete the key entirely, as __getstate__ must return a non-empty
        # dict in order for our __setstate__ method to be called).
        d = self.__dict__.copy()
        for key in self.__dict__.keys():
            if key.endswith("_modpt"):
                d[key] = None
        return d


# Modeller 9 compatibility
class modobject(ModObject):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(modobject)
        ModObject.__init__(self, *args, **keys)
