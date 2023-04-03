# Example for: Selection.energy()

# This will calculate the stereochemical energy (bonds,
# angles, dihedrals, impropers) for a given model.

from modeller import *
from modeller.scripts import complete_pdb

env = Environ()
env.io.atom_files_directory = ['../atom_files']
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

def patch_disulfides(mdl):
    # Must patch disulfides here to calculate the non-bonded
    # energy properly. Also, when you use hydrogens, disulfides
    # must always be patched so that sulfhydril hydrogens are
    # removed from the model.
    for ids in [ ('17:A', '39:A'),
                 ( '3:A', '22:A'),
                 ('53:A', '59:A'),
                 ('41:A', '52:A') ]:
        mdl.patch(residue_type='DISU', residues=[mdl.residues[r] for r in ids])

mdl = complete_pdb(env, "1fas", special_patches=patch_disulfides)

# Select all atoms
atmsel = Selection(mdl)

mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)

# Actually calculate the energy
(molpdf, terms) = atmsel.energy(edat=EnergyData(dynamic_sphere=True))

# molpdf is the total 'energy', and terms contains the contributions from
# each physical type. Here we print out the bond length contribution:
print("Bond energy is %.3f" % terms[physical.bond])
