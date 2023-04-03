# Illustrates the SALIGN multiple sequence alignment
from modeller import *

log.verbose()
env = Environ()
env.io.atom_files_directory = ['.', '../atom_files']

aln = Alignment(env, file='malign_in.ali')

aln.salign(overhang=30, gap_penalties_1d=(-450, -50),
           alignment_type='tree', output='ALIGNMENT')

aln.write(file='malign.ali', alignment_format='PIR')
