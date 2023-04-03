"""Classes for loop modeling using the DOPE potential"""

from modeller.automodel.loopmodel import LoopModel
from modeller import GroupRestraints, physical, gbsa
from modeller.util.deprecation import _deprecation_handler


__docformat__ = "epytext en"


class DOPELoopModel(LoopModel):
    """Loop modeling using the DOPE potential"""

    def __init__(self, env, sequence, alnfile=None, knowns=None, inimodel=None,
                 deviation=None, library_schedule=None, csrfile=None,
                 inifile=None, assess_methods=None, loop_assess_methods=None,
                 root_name=None):
        LoopModel.__init__(self, env, sequence, alnfile, knowns, inimodel,
                           deviation, library_schedule, csrfile, inifile,
                           assess_methods, loop_assess_methods, root_name)
        self.loop.env.schedule_scale = physical.Values(default=1.0,
                                                       nonbond_spline=0.6)
        edat = self.loop.env.edat
        edat.contact_shell = 8.00
        edat.dynamic_sphere = False
        edat.dynamic_lennard = True
        edat.dynamic_coulomb = False
        edat.relative_dielectric = 1.0
        edat.dynamic_modeller = True
        edat.energy_terms.append(gbsa.Scorer(cutoff=edat.contact_shell))

    def read_potential(self):
        return GroupRestraints(self.env, classes='$(LIB)/atmcls-mf.lib',
                               parameters='$(LIB)/dist-mf.lib')

    def loop_restraints(self, atmsel, aln):
        dih_lib_only = True
        mnch_lib = 1
        res_at = 1
        self.restraints.clear()
        for typ in ('bond', 'angle', 'improper', 'dihedral'):
            self.restraints.make(atmsel, aln=aln, restraint_type=typ,
                                 spline_on_site=self.spline_on_site,
                                 dih_lib_only=dih_lib_only,
                                 mnch_lib=mnch_lib, restraint_sel_atoms=res_at)
        for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
            self.restraints.make(atmsel, aln=aln,
                                 restraint_type=typ+'_dihedral',
                                 spline_on_site=self.spline_on_site,
                                 dih_lib_only=dih_lib_only, mnch_lib=mnch_lib,
                                 restraint_sel_atoms=res_at, spline_range=4.0,
                                 spline_dx=0.3, spline_min_points=5)
        # Generate restraints on all non-standard residues:
        self.nonstd_restraints(aln)

        self.special_restraints(aln)


# Modeller 9 compatibility
class dope_loopmodel(DOPELoopModel):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(dope_loopmodel)
        DOPELoopModel.__init__(self, *args, **keys)
