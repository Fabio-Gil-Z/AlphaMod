# Example for: Model.orient()

# This will orient the model along the principal axes of the inertia ellipsoid:

from modeller import *

env = Environ()
env.io.atom_files_directory = ['../atom_files']
mdl = Model(env)
mdl.read(file='1fas')
r = mdl.orient()
mdl.write(file='1fas.ini')

print("Translation: " + str(r.translation))
