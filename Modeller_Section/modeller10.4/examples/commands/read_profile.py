# Example file for: Profile.read(), Profile.to_alignment()

from modeller import *

env = Environ()

# Create a new, blank, profile
prf = Profile(env)

# Read in the profile file
prf.read(file='toxin.prf', profile_format='TEXT')

# Convert the profile to alignment
aln = prf.to_alignment()

# Write out the alignment
aln.write(file='readprofile.pir', alignment_format='PIR')
