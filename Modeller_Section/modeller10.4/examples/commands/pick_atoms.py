# This will pick various subsets of atoms in the MODEL and compare them
# with MODEL2.

from modeller import *

env = Environ()
env.io.atom_files_directory = ['../atom_files']
log.level(1, 1, 1, 1, 0)

# Read the models and the alignment:
mdl  = Model(env, file='1fas')
mdl2 = Model(env, file='2ctx')
aln = Alignment(env, file='toxin.ali', align_codes=('1fas', '2ctx'))
aln.write(file='toxin.pap', alignment_format='PAP')

# Pick and superpose mainchain atoms:
atmsel = Selection(mdl).only_mainchain()
atmsel.superpose(mdl2, aln)

# Pick and superpose sidechain atoms:
atmsel = Selection(mdl).only_sidechain()
atmsel.superpose(mdl2, aln)

# Pick and superpose CA and CB atoms:
atmsel = Selection(mdl).only_atom_types('CA CB')
atmsel.superpose(mdl2, aln)

# Pick all atoms in residues with rings:
atmsel = Selection(mdl).only_residue_types('TYR PHE TRP HIS')

# Pick and superpose all atoms:
atmsel = Selection(mdl)
atmsel.superpose(mdl2, aln)

# Pick and superpose CA and CB atoms in one segment only:
atmsel = Selection(mdl.residue_range('2:A', '10:A')).only_atom_types('CA CB')
atmsel.superpose(mdl2, aln)

# Pick and superpose all atoms within 6 angstroms of the 'CA' atom in
# residue '10' in chain A:
atmsel = mdl.atoms['CA:10:A'].select_sphere(6.0)
atmsel.superpose(mdl2, aln)

# Pick and superpose all atoms within 6 angstroms of any atom in
# segment 2:A to 10:A
atmsel = Selection(mdl.residue_range('2:A', '10:A')).select_sphere(6.0)
atmsel.superpose(mdl2, aln)

# Pick all atoms in the model
atmsel = Selection(mdl)

# Pick all atoms in all loops (ie residues within 2 positions
# of any gap in the alignment):
loops = mdl2.loops(aln, minlength=5, maxlength=15, insertion_ext=2,
                   deletion_ext=2)
atmsel = Selection(loops)

# Pick all atoms within 6 angstroms of all loops
atmsel = Selection(loops).select_sphere(6.0)
