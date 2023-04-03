# Comparative modeling by the AutoModel class, generating mmCIF outputs
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the AutoModel class

log.verbose()    # request verbose output
env = Environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = AutoModel(env,
              alnfile  = 'alignment.ali',     # alignment filename
              knowns   = '5fd1',              # codes of the templates
              sequence = '1fdx')              # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
                                    # (determines how many models to calculate)
a.set_output_model_format("MMCIF")  # request mmCIF rather than PDB outputs
a.make()                            # do the actual comparative modeling
