# Example for: Model.res_num_from()

# This will transfer residue numbers and chain ids from model2 to model.

from modeller import *

log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
env = Environ()
env.io.atom_files_directory = ['../atom_files']

# Read an alignment for the transfer
aln = Alignment(env, file='toxin.ali', align_codes=('2ctx', '1fas'))
# Read the template and target models:
mdl2 = Model(env, file='2ctx')
mdl  = Model(env, file='1fas')
# Transfer the residue and chain ids and write out the new MODEL:
mdl.res_num_from(mdl2, aln)
mdl.write(file='1fas.ini')
