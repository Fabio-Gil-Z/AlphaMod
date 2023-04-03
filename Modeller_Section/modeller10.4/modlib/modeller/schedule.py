"""Classes for defining an optimization schedule"""

from modeller import physical
from modeller.util.deprecation import _deprecation_handler

__docformat__ = "epytext en"


class Schedule(list):
    """A set of optimization steps (schedule)"""

    def __init__(self, last_scales, steps):
        self.last_scales = last_scales
        list.__init__(self, steps)

    def make_for_model(self, mdl):
        """Returns a new schedule, trimmed to the right length for the model"""
        reslen = len(mdl.residues)
        obj = Schedule(self.last_scales, ())
        for x in self:
            obj.append(x.copy())
            if x.max_residue is not None and x.max_residue != 9999 \
               and x.max_residue >= reslen:
                break
        for x in range(0, self.last_scales):
            if x >= len(obj) or x >= len(self):
                break
            obj[-1-x].scale = self[-1-x].scale
        return obj

    def __mul__(self, other):
        obj = Schedule(self.last_scales, ())
        for x in self:
            obj.append(x * other)
        return obj

    def write(self, fh):
        """Write out a schedule to a file, in legacy format"""
        nphycns = len(self[0].scale)
        fh.write("# MTH NRNG"
                 + "".join(["%6d" % (x+1) for x in range(nphycns)]) + "\n")
        fh.write("-" * (11+nphycns*6) + "\n")
        for (num, step) in enumerate(self):
            mr = step.max_residue
            if mr is None:
                mr = 9999
            fh.write("%2d%3d%5d " % (num+1, 1, mr)
                     + "".join(["%6.3f" % x for x in step.scale]) + "\n")


class Step(object):
    """A single step in an optimization schedule"""

    def __init__(self, optimizer, max_residue, scale):
        self.optimizer = optimizer
        self.max_residue = max_residue
        self.scale = scale

    def copy(self):
        return Step(self.optimizer, self.max_residue, self.scale)

    def __mul__(self, other):
        if not isinstance(other, physical.Values):
            raise TypeError("Must use physical.Values objects here")
        return Step(self.optimizer, self.max_residue, self.scale * other)

    def optimize(self, atmsel, **vars):
        mdl = atmsel.get_model()
        if self.max_residue is not None:
            span = (0, self.max_residue)
            mdl.restraints.unpick_all()
            mdl.restraints.pick(atmsel, residue_span_range=span)
            vars['residue_span_range'] = span
        o = self.optimizer
        if isinstance(o, int):
            raise TypeError(
                "Should use real optimizer objects now, not 1 or 3")
        if o.__class__ == type:
            o = o()   # create an instance if a class object was used
        return o.optimize(atmsel, schedule_scale=self.scale, **vars)


# Modeller 9 compatibility
class schedule(Schedule):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(schedule)
        Schedule.__init__(self, *args, **keys)


class step(Step):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(step)
        Step.__init__(self, *args, **keys)
