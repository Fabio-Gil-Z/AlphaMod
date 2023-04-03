# Example for: Selection.superpose()

# This will use a given alignment to superpose Calpha atoms of
# one structure (2ctx) on the other (1fas).

from modeller import *

env = Environ()
env.io.atom_files_directory = ['../atom_files']

mdl  = Model(env, file='1fas')
mdl2 = Model(env, file='2ctx')
aln = Alignment(env, file='toxin.ali', align_codes=('1fas', '2ctx'))

atmsel = Selection(mdl).only_atom_types('CA')
r = atmsel.superpose(mdl2, aln)

# We can now use the calculated RMS, DRMS, etc. from the returned 'r' object:
rms = r.rms
drms = r.drms
print("%d equivalent positions" % r.num_equiv_pos)

mdl2.write(file='2ctx.fit')
