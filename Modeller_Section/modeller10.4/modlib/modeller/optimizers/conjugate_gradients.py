import _modeller
from modeller.optimizers.builtin_optimizer import BuiltinOptimizer
from modeller.util.deprecation import _deprecation_handler


class ConjugateGradients(BuiltinOptimizer):
    """Optimize with Beale restraint conjugate gradients method"""
    _modpt = None
    __free_func = _modeller.mod_cg_optimizer_free
    __ok_keys = ('edat', 'libs', 'schedule_scale', 'residue_span_range',
                 'max_iterations', 'min_atom_shift', 'output', 'actions')

    def __init__(self, output='NO_REPORT', min_atom_shift=0.01,
                 residue_span_range=(0, 99999), **vars):
        BuiltinOptimizer.__init__(self)
        self._modpt = _modeller.mod_cg_optimizer_new(self)
        self.__params = {}
        for key in ("output", "min_atom_shift", "residue_span_range"):
            vars[key] = eval(key)
        self._update_params(self.__params, self.__ok_keys, vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._modpt = _modeller.mod_cg_optimizer_new(self)

    def __del__(self):
        if self._modpt and _modeller:
            self.__free_func(self._modpt)

    def __get_optpt(self):
        return _modeller.mod_cg_optimizer_opt_get(self._modpt)

    def optimize(self, atmsel, **vars):
        """Optimize the given atom selection"""
        (mdl, libs, edat, inds, vars) = \
            self._prep_builtin_optimizer(atmsel, self.__params, self.__ok_keys,
                                         vars)
        self._prep_actions(vars)
        self.atmsel = atmsel
        func = _modeller.mod_cg_optimize
        ret = func(self._modpt, mdl.modpt, edat.modpt, libs.modpt,
                   inds, **vars)
        self.atmsel = None
        return ret

    def trace(self, out):
        """Write out information about the current optimization step"""
        if self.step == 0:
            out.write("# Conjugate gradients optimization\n")
            out.write("# %5s %16s %8s %8s %7s %16s\n"
                      % ('Step', 'Current energy', 'Av shift',
                         'Mx shift', 'Funcs', 'Gradient'))
        log = "%7d %16.5f %8.4f %8.4f %7d %16.6g\n" \
              % (self.step, self.current_e, self.shiftavr,
                 self.shiftmax, self.funcs, self.gsq)
        out.write(log)

    def __get_gsq(self):
        return _modeller.mod_cg_optimizer_gsq_get(self._modpt)

    optpt = property(__get_optpt)
    gsq = property(__get_gsq, doc="Current gradient")


# Modeller 9 compatibility
class conjugate_gradients(ConjugateGradients):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(conjugate_gradients)
        ConjugateGradients.__init__(self, *args, **keys)
