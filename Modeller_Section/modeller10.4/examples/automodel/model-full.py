# A sample script for fully automated comparative modeling
from modeller import *
from modeller.automodel import *    # Load the AutoModel class

log.verbose()
env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = AutoModel(env,
              # file with template codes and target sequence
              alnfile  = 'alignment.seg',
              # PDB codes of the templates
              knowns   = ('5fd1', '1fdn', '1fxd', '1iqz'),
              # code of the target
              sequence = '1fdx')
a.auto_align()                      # get an automatic alignment
a.make()                            # do comparative modeling
