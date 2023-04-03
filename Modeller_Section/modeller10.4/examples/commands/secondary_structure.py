# Example for Model.build_sequence(), secondary_structure.Alpha()

from modeller import *
from modeller.optimizers import ConjugateGradients

# Set up environment
e = Environ()
e.libs.topology.read('${LIB}/top_heav.lib')
e.libs.parameters.read('${LIB}/par.lib')

# Build an extended chain model from primary sequence, and write it out
m = Model(e)
m.build_sequence('GSCASVCGV')
m.write(file='extended-chain.pdb')

# Make stereochemical restraints on all atoms
allatoms = Selection(m)
m.restraints.make(allatoms, restraint_type='STEREO', spline_on_site=False)

# Constrain all residues to be alpha-helical
# (Could also use m.residue_range() rather than m.residues here.)
m.restraints.add(secondary_structure.Alpha(m.residues))

# Get an optimized structure with CG, and write it out
cg = ConjugateGradients()
cg.optimize(allatoms, max_iterations=100)
m.write(file='alpha-helix.pdb')
