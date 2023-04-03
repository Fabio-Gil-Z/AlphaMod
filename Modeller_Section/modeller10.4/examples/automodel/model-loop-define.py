from modeller import *
from modeller.automodel import *

log.verbose()
env = Environ()

env.io.atom_files_directory = ['.', '../atom_files']

# Create a new class based on 'LoopModel' so that we can redefine
# select_loop_atoms
class MyLoop(LoopModel):
    # This routine picks the residues to be refined by loop modeling
    def select_loop_atoms(self):
        # Two residue ranges (both will be refined simultaneously)
        return Selection(self.residue_range('19:A', '28:A'),
                         self.residue_range('45:A', '50:A'))

a = MyLoop(env,
           alnfile  = 'alignment.ali',      # alignment filename
           knowns   = '5fd1',               # codes of the templates
           sequence = '1fdx',               # code of the target
           loop_assess_methods=assess.DOPE) # assess each loop with DOPE
a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 2           # Last loop model

a.make()                            # do modeling and loop refinement
