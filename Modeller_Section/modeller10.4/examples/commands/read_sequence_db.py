from modeller import *

log.verbose()
env = Environ()

sdb = SequenceDB(env)
sdb.read(seq_database_file='pdb95.fsa', seq_database_format='FASTA',
         chains_list='ALL', minmax_db_seq_len=[1, 40000], clean_sequences=True)

sdb.write(seq_database_file='pdb95.bin', seq_database_format='BINARY',
          chains_list='ALL')
