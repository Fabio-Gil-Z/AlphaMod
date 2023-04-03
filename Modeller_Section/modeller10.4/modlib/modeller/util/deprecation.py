import _modeller
from modeller.util.logger import log


class _DeprecationHandler(object):
    def __init__(self):
        self.depmode = _modeller.mod_deprecation_mode_get()

    error = property(lambda self: self.depmode == 'ERROR')

    def _init_class(self, cls):
        clsname = cls.__name__
        basename = cls.__bases__[0].__name__
        msg = ("The class '%s' is deprecated; use '%s' instead"
               % (clsname, basename))
        if self.error:
            log.error(clsname, msg)
        else:
            log.warning(clsname, msg)


_deprecation_handler = _DeprecationHandler()
