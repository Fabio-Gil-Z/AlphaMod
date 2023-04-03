# Example for PSF and binary trajectory output

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import MolecularDynamics, actions

env = Environ()
env.io.atom_files_directory = ['../atom_files']
env.edat.dynamic_sphere = True
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

code = '1fas'
mdl = complete_pdb(env, code)

# Stereochemical restraints on all atoms:
atmsel = Selection(mdl)
mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)

# Write a PSF
mdl.write_psf(code+'.psf')

# Run 100 steps of MD, writing a CHARMM binary trajectory every 5 steps
md = MolecularDynamics(output='REPORT')
md.optimize(atmsel, temperature=300, max_iterations=100,
            actions=actions.CHARMMTrajectory(5, filename=code+'.dcd'))
