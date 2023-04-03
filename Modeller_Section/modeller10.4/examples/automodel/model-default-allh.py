from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()

env.io.atom_files_directory = ['.', '../atom_files']

a = AllHModel(env, alnfile='alignment.ali', knowns='5fd1', sequence='1fdx')
a.starting_model = a.ending_model = 4

a.make()
