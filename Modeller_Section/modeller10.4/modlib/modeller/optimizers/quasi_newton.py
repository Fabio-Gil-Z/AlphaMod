import _modeller
from modeller.optimizers.builtin_optimizer import BuiltinOptimizer
from modeller.util.deprecation import _deprecation_handler


class QuasiNewton(BuiltinOptimizer):
    """Optimize with Broyden-Fletcher-Goldfarb-Shanno variable metric
       (quasi Newton) method"""
    _modpt = None
    __free_func = _modeller.mod_qn_optimizer_free
    __ok_keys = ('edat', 'libs', 'schedule_scale', 'residue_span_range',
                 'max_iterations', 'min_atom_shift', 'max_atom_shift',
                 'output', 'actions')

    def __init__(self, output='NO_REPORT', min_atom_shift=0.01,
                 max_atom_shift=100.0, residue_span_range=(0, 99999), **vars):
        BuiltinOptimizer.__init__(self)
        self._modpt = _modeller.mod_qn_optimizer_new(self)
        self.__params = {}
        for key in ("output", "min_atom_shift", "max_atom_shift",
                    "residue_span_range"):
            vars[key] = eval(key)
        self._update_params(self.__params, self.__ok_keys, vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._modpt = _modeller.mod_qn_optimizer_new(self)

    def __del__(self):
        if self._modpt:
            self.__free_func(self._modpt)

    def __get_optpt(self):
        return _modeller.mod_qn_optimizer_opt_get(self._modpt)

    def optimize(self, atmsel, **vars):
        """Optimize the given atom selection"""
        (mdl, libs, edat, inds, vars) = \
            self._prep_builtin_optimizer(atmsel, self.__params, self.__ok_keys,
                                         vars)
        self._prep_actions(vars)
        self.atmsel = atmsel
        func = _modeller.mod_qn_optimize
        ret = func(self._modpt, mdl.modpt, edat.modpt, libs.modpt,
                   inds, **vars)
        self.atmsel = None
        return ret

    def trace(self, out):
        """Write out information about the current optimization step"""
        if self.step == 0:
            out.write("# Quasi Newton optimization\n")
            out.write("# %5s %16s %8s %8s %7s\n"
                      % ('Step', 'Current energy', 'Av shift',
                         'Mx shift', 'Funcs'))
        log = "%7d %16.5f %8.4f %8.4f %7d\n" \
              % (self.step, self.current_e, self.shiftavr,
                 self.shiftmax, self.funcs)
        out.write(log)

    optpt = property(__get_optpt)


# Modeller 9 compatibility
class quasi_newton(QuasiNewton):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(quasi_newton)
        QuasiNewton.__init__(self, *args, **keys)
