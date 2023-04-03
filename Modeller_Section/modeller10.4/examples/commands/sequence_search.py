# Example for: SequenceDB.search()

# This will search the MODELLER database of representative protein chains
# for chains similar to the specified sequence.

from modeller import *

log.verbose()
env = Environ()

# Read in the sequences of all PDB structures
try:
    sdb = SequenceDB(env, seq_database_file='pdball.pir',
                     seq_database_format='PIR',
                     chains_list='very-short-for-test.cod')
except IOError:
    print("""
Could not read sequence database file. This file is not included by default
in the Modeller distribution, but you can download it from the Modeller
downloads page (https://salilab.org/modeller/supplemental.html).

Note: it is recommended to use Profile.build() rather than SequenceDB.search().
See step 1 of the Modeller basic tutorial at
https://salilab.org/modeller/tutorial/basic.html
""")
    raise

# Read in the query sequence in alignment format
aln = Alignment(env, file='toxin.ali', align_codes='2nbt')

sdb.search(aln, search_randomizations=20, # should use 100 in real life
           seq_database_file='pdball.pir',
           search_group_list='pdb_95.grp',
           off_diagonal=9999, gap_penalties_1d=(-800, -400),
           signif_cutoff=(1.5, 5.0))

aln.malign()
aln.write(file='toxin-search.pap', alignment_format='PAP')
