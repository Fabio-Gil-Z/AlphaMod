# Comparative modeling by the AutoModel class
from modeller import *
from modeller.automodel import *    # Load the AutoModel class

log.verbose()
env = Environ()

class NewLoop(LoopModel):
    def select_loop_atoms(self):
        return Selection(self.residue_range('2:A', '5:A'))


# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = NewLoop(env,
            alnfile  = 'alignment-tiny.ali',     # alignment filename
            knowns   = '5fd1',                   # codes of the templates
            sequence = '1fdx')                   # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 3                 # index of the last model
                                    # (determines how many models to calculate)
a.md_level = refine.very_fast
a.max_var_iterations = 50
a.final_malign3d = True

a.loop.starting_model = 1
a.loop.ending_model   = 2
a.loop.md_level       = refine.very_fast

a.make()                            # do comparative modeling
