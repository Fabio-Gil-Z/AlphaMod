# Comparative modeling by the AutoModel class
#
# Demonstrates how to refine only a part of the model.
#
# You may want to use the more exhaustive "loop" modeling routines instead.
#
from modeller import *
from modeller.automodel import *    # Load the AutoModel class

log.verbose()

# Override the 'select_atoms' routine in the 'AutoModel' class:
# (To build an all-hydrogen model, derive from AllHModel rather than AutoModel
# here.)
class MyModel(AutoModel):
    def select_atoms(self):
        # Select residues 1 and 2 in chain A (PDB numbering)
        return Selection(self.residue_range('1:A', '2:A'))

        # Residues 4, 6, 10 in chain A:
        # return Selection(self.residues['4:A'], self.residues['6:A'],
        #                  self.residues['10:A'])

        # All residues except 1-5 in chain A:
        # return Selection(self) - Selection(self.residue_range('1:A', '5:A'))


env = Environ()
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
# selected atoms do not feel the neighborhood
env.edat.nonbonded_sel_atoms = 2

# Be sure to use 'MyModel' rather than 'AutoModel' here!
a = MyModel(env,
            alnfile  = 'alignment.ali',     # alignment filename
            knowns   = '5fd1',              # codes of the templates
            sequence = '1fdx')              # code of the target

a.starting_model= 3                # index of the first model
a.ending_model  = 3                # index of the last model
                                   # (determines how many models to calculate)
a.make()                           # do comparative modeling
