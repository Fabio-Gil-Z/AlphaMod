"""Optimization schedules used by AutoModel and LoopModel in the initial
   conjugate gradients optimization (VTFM).

   Note: if adapting these schedules for your own uses, it is recommended
   that the scaling factors for the final step are always set to
   physical.Values(default=1.00). This is because AutoModel scales energies
   by both these factors and schedule_scale during VTFM, but only by
   schedule_scale during MD&SA refinement."""

from modeller import physical
from modeller.schedule import Schedule, Step
from modeller.optimizers import ConjugateGradients as CG

__docformat__ = "epytext en"


def mk_scale(default, nonbond, spline=None):
    """Utility function for generating scaling values"""
    v = physical.Values(default=default)
    for term in (physical.soft_sphere, physical.lennard_jones,
                 physical.coulomb, physical.gbsa, physical.em_density,
                 physical.saxs):
        v[term] = nonbond
    if spline is not None:
        v[physical.nonbond_spline] = spline
    return v


#: thorough optimization
slow = Schedule(
    4, [Step(CG, 2, mk_scale(default=0.01, nonbond=0.0)),
        Step(CG, 4, mk_scale(default=0.10, nonbond=0.0)),
        Step(CG, 6, mk_scale(default=0.50, nonbond=0.0))] +
       [Step(CG, rng, mk_scale(default=1.00, nonbond=0.0)) for rng in
        (8, 10, 14, 18, 20, 24, 30, 25, 40, 45, 50, 55, 60, 70, 80, 90, 100,
         120, 140, 160, 200, 250, 300, 400, 500)] +
       [Step(CG, 600, mk_scale(default=1.00, nonbond=0.01)),
        Step(CG, 800, mk_scale(default=1.00, nonbond=0.1)),
        Step(CG, 1000, mk_scale(default=1.00, nonbond=0.5)),
        Step(CG, 9999, physical.Values(default=1.00))])

#: normal optimization
normal = Schedule(
    4,
    [Step(CG, 2, mk_scale(default=0.01, nonbond=0.0)),
     Step(CG, 4, mk_scale(default=0.10, nonbond=0.0)),
     Step(CG, 6, mk_scale(default=0.50, nonbond=0.0))] +
    [Step(CG, rng, mk_scale(default=1.00, nonbond=0.0)) for rng in
     (10, 20, 30, 50, 80, 120, 200, 300)] +
    [Step(CG, 500, mk_scale(default=1.00, nonbond=0.01)),
     Step(CG, 800, mk_scale(default=1.00, nonbond=0.1)),
     Step(CG, 1000, mk_scale(default=1.00, nonbond=0.5)),
     Step(CG, 9999, physical.Values(default=1.00))])

#: fast optimization
fast = Schedule(
    4, [Step(CG, 2, mk_scale(default=0.01, nonbond=0.0)),
        Step(CG, 6, mk_scale(default=0.10, nonbond=0.0)),
        Step(CG, 10, mk_scale(default=0.50, nonbond=0.0))] +
       [Step(CG, rng, mk_scale(default=1.00, nonbond=0.0)) for rng in
        (20, 50, 100, 200)] +
       [Step(CG, 500, mk_scale(default=1.00, nonbond=0.01)),
        Step(CG, 800, mk_scale(default=1.00, nonbond=0.1)),
        Step(CG, 1000, mk_scale(default=1.00, nonbond=0.5)),
        Step(CG, 9999, physical.Values(default=1.00))])

#: very fast optimization
very_fast = Schedule(
    0, [Step(CG, 9999, mk_scale(default=0.01, nonbond=0.0)),
        Step(CG, 9999, mk_scale(default=0.10, nonbond=0.0)),
        Step(CG, 9999, mk_scale(default=0.50, nonbond=0.0)),
        Step(CG, 9999, physical.Values(default=1.00, soft_sphere=0.01,
                                       lennard_jones=0, coulomb=0, gbsa=0.01,
                                       em_density=0.01, saxs=0.01)),
        Step(CG, 9999, mk_scale(default=1.00, nonbond=0.1)),
        Step(CG, 9999, mk_scale(default=1.00, nonbond=0.5)),
        Step(CG, 9999, physical.Values(default=1.00))])

#: fastest optimization
fastest = Schedule(
    3, [Step(CG, 2, mk_scale(default=0.01, nonbond=0.0)),
        Step(CG, 5, mk_scale(default=0.50, nonbond=0.0)),
        Step(CG, 10, mk_scale(default=1.00, nonbond=0.0)),
        Step(CG, 50, mk_scale(default=1.00, nonbond=0.0)),
        Step(CG, 200, mk_scale(default=1.00, nonbond=0.0)),
        Step(CG, 600, mk_scale(default=1.00, nonbond=0.01)),
        Step(CG, 1000, mk_scale(default=1.00, nonbond=0.50)),
        Step(CG, 9999, physical.Values(default=1.00))])

#: optimization with no residue span restriction (good for loop modeling)
loop = Schedule(
    4, [Step(CG, None, mk_scale(default=1.00, nonbond=0.0, spline=1.00)),
        Step(CG, None, mk_scale(default=1.00, nonbond=0.01, spline=0.01)),
        Step(CG, None, mk_scale(default=1.00, nonbond=0.10, spline=0.10)),
        Step(CG, None, mk_scale(default=1.00, nonbond=0.50, spline=0.50)),
        Step(CG, None, physical.Values(default=1.00))])
