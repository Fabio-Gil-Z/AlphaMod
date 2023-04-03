# This will swap certain atom names in some planar sidechains to satisfy
# the IUPAC convention.

from modeller import *

env = Environ()
env.io.atom_files_directory = ['../atom_files']
log.level(1, 1, 1, 1, 0)

mdl = Model(env, file='2abx')
mdl.to_iupac()
mdl.write(file='2abx.iup')
