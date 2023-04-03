from modeller.library_restraints import bonds, angles, impropers
from modeller.library_restraints import omega_dihedrals, phi_psi
from modeller.library_restraints import chi1, chi2, chi3, chi4

def make_restraints(atmsel, restraints, num_selected):
    for a in (bonds, angles, impropers, omega_dihedrals, phi_psi,
              chi1, chi2, chi3, chi4):
        a.make_restraints(atmsel, restraints, num_selected)
