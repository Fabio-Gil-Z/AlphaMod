# Example for: Alignment.malign()

# This will read all sequences from a file, align them, and write
# the alignment to a new file:

from modeller import *

env = Environ()

aln = Alignment(env, file='toxin.ali', align_codes='all')
aln.malign(gap_penalties_1d=(-600, -400))
aln.write(file='toxin-seq.pap', alignment_format='PAP')
