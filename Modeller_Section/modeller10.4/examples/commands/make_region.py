# Example for: Model.make_region()

# This will define a random contiguous patch of atoms on a surface of the
# protein.

from modeller import *

env = Environ(rand_seed=-18343)
log.level(1, 1, 1, 1, 0)

# Read the PDB file
mdl = Model(env)
mdl.read(file='../atom_files/pdb1fdn.ent')

# Calculate atomic accessibilities (in Biso) with appropriate probe_radius
myedat = EnergyData()
myedat.radii_factor = 1.6
mdl.write_data(edat=myedat, output='PSA ATOMIC_SOL',
               psa_integration_step=0.05, probe_radius=0.1)

# Get the "random" patch of exposed atoms on the surface
mdl.make_region(atom_accessibility=0.5, region_size=35)

# Write out a PDB file with the patch indicated by Biso = 1:
mdl.write(file='1fdn.reg')

# Can also select the patch residues and use selection methods:
s = Selection([a for a in mdl.atoms if a.biso > 0.99])
print("%d atoms in surface patch" % len(s))
