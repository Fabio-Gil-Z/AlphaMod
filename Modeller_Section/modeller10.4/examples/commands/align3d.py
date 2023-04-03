# Example for: Alignment.align3d(), Selection.superpose()

# This will align 3D structures of two proteins:

from modeller import *
log.verbose()
env = Environ()
env.io.atom_files_directory = ['../atom_files']

# First example: read sequences from a sequence file:
aln = Alignment(env)
aln.append(file='toxin.ali', align_codes=['1fas', '2ctx'])
aln.align(gap_penalties_1d=[-600, -400])
aln.align3d(gap_penalties_3d=[0, 4.0])
aln.write(file='toxin-str.ali')

# Second example: read sequences from PDB files to eliminate the
# need for the toxin.ali sequence file:
mdl = Model(env)
aln = Alignment(env)
for code in ['1fas', '2ctx']:
    mdl.read(file=code)
    aln.append_model(mdl, align_codes=code, atom_files=code)
aln.align(gap_penalties_1d=(-600, -400))
aln.align3d(gap_penalties_3d=(0, 2.0))
aln.write(file='toxin-str.ali')

# And now superpose the two structures using current alignment to get
# various RMS's:
mdl = Model(env, file='1fas')
atmsel = Selection(mdl).only_atom_types('CA')
mdl2 = Model(env, file='2ctx')
atmsel.superpose(mdl2, aln)
