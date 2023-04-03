# This will create a VTFM optimization schedule and then
# use it to optimize the model

from modeller import *
from modeller.scripts import complete_pdb

# Load in optimizer and schedule support
from modeller import schedule, optimizers

log.verbose()

env = Environ()
env.io.atom_files_directory = ['../atom_files']
env.edat.dynamic_sphere = True

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')
code = '1fas'
mdl = complete_pdb(env, code)

# Generate the restraints:
atmsel = Selection(mdl)
mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)

# Create our own library schedule:
# 5 steps of conjugate gradients (CG), each step using a larger
# residue range (2 up to 9999) and energy scaling factor (0.01 up to 1.0),
# followed by 3 steps of molecular dynamics (MD) at successively lower
# temperature. The scaling factors for the last 5 steps are always retained.
CG = optimizers.ConjugateGradients
MD = optimizers.MolecularDynamics
libsched = schedule.Schedule(5,
          [ schedule.Step(CG, 2, physical.Values(default=0.01)),
            schedule.Step(CG, 5, physical.Values(default=0.1)),
            schedule.Step(CG, 10, physical.Values(default=0.2)),
            schedule.Step(CG, 50, physical.Values(default=0.5)),
            schedule.Step(CG, 9999, physical.Values(default=1.0)),
            schedule.Step(MD(temperature=300.), 9999, \
                          physical.Values(default=1.0)),
            schedule.Step(MD(temperature=200.), 9999, \
                          physical.Values(default=1.0)),
            schedule.Step(MD(temperature=100.), 9999, \
                          physical.Values(default=1.0)) ])

# Make a trimmed schedule suitable for our model, and scale it by schedule_scale
mysched = libsched.make_for_model(mdl) * env.schedule_scale

# Write the trimmed schedule to a file
fh = open(code+'.sch', 'w')
mysched.write(fh)
fh.close()

# Optimize for all steps in the schedule
for step in mysched:
    step.optimize(atmsel, output='REPORT', max_iterations=200)

mdl.write(file=code+'.B')
