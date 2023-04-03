# Example for: Selection.debug_function()

# This will use the MODELLER AutoModel class to construct homology
# restraints for 1fas. It will then use Model.debug_function() to test
# the source code for the function and derivatives calculation
# by comparing analytical and numerical first derivatives (note that
# AutoModel is a derived class of model, so all 'model' methods will work
# on 'AutoModel'). Some discrepancies may be reported but ignore them here.

from modeller import *
from modeller.automodel import AutoModel    # Load the AutoModel class

log.verbose()
env = Environ()
env.io.atom_files_directory = ['../atom_files']

a = AutoModel(env, alnfile = 'debug_function.ali',
              knowns  = ('2ctx', '2nbt'), sequence = '1fas')
a.spline_on_site = False
a.make(exit_stage=1)

# Test on all atoms
atmsel = Selection(a)

# To assign 0 weights to restraints whose numerical derivatives
# code does not work (i.e., splines for angles and dihedrals):
scal = physical.Values(default=1.0, lennard_jones=0, coulomb=0, h_bond=0,
                       phi_dihedral=0, psi_dihedral=0, omega_dihedral=0,
                       chi1_dihedral=0, chi2_dihedral=0, chi3_dihedral=0,
                       chi4_dihedral=0, disulfide_angle=0,
                       disulfide_dihedral=0, chi5_dihedral=0)
atmsel.energy(output='SHORT', schedule_scale=scal)
atmsel.debug_function(debug_function_cutoff=(15.00, 0.10, 0.1),
                      detailed_debugging=True, schedule_scale=scal)
