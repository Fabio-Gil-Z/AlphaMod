from modeller.util.deprecation import _deprecation_handler


class Optimizer(object):
    atmsel = None

    def optimize(self, atmsel):
        raise NotImplementedError

    def get_selection(self):
        return self.atmsel

    def _update_params(self, params, ok_keys, vars):
        for key in vars.keys():
            if key in ok_keys:
                params[key] = vars[key]
            else:
                raise ValueError("Unrecognized parameter: %s" % key)


# Modeller 9 compatibility
class optimizer(Optimizer):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(optimizer)
        Optimizer.__init__(self, *args, **keys)
