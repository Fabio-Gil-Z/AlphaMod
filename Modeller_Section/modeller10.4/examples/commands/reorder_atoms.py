# Example for: Model.reorder_atoms()

# This will standardize the order of atoms in the model.

from modeller import *

env = Environ()
env.io.atom_files_directory = ['../atom_files']

# Order the atoms according to a topology library:
env.libs.topology.read(file='$(LIB)/top_heav.lib')

mdl = Model(env, file='1fas')
mdl.reorder_atoms()
mdl.write(file='1fas.ini1')
