# Example for: ConjugateGradients(), MolecularDynamics(), Model.switch_trace()

# This will optimize stereochemistry of a given model, including
# non-bonded contacts.

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import ConjugateGradients, MolecularDynamics, actions

env = Environ()
env.io.atom_files_directory = ['../atom_files']
env.edat.dynamic_sphere = True

env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

code = '1fas'
mdl = complete_pdb(env, code)
mdl.write(file=code+'.ini')

# Select all atoms:
atmsel = Selection(mdl)

# Generate the restraints:
mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)
mdl.restraints.write(file=code+'.rsr')

mpdf = atmsel.energy()

# Create optimizer objects and set defaults for all further optimizations
cg = ConjugateGradients(output='REPORT')
md = MolecularDynamics(output='REPORT')

# Open a file to get basic stats on each optimization
trcfil = open(code+'.D00000001', 'w')

# Run CG on the all-atom selection; write stats every 5 steps
cg.optimize(atmsel, max_iterations=20, actions=actions.Trace(5, trcfil))
# Run MD; write out a PDB structure (called '1fas.D9999xxxx.pdb') every
# 10 steps during the run, and write stats every 10 steps
md.optimize(atmsel, temperature=300, max_iterations=50,
            actions=[actions.WriteStructure(10, code+'.D9999%04d.pdb'),
                     actions.Trace(10, trcfil)])
# Finish off with some more CG, and write stats every 5 steps
cg.optimize(atmsel, max_iterations=20,
            actions=[actions.Trace(5, trcfil)])

mpdf = atmsel.energy()

mdl.write(file=code+'.B')
