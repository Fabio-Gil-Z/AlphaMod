from modeller import *
from modeller.scripts import complete_pdb

# Load the C extension module; this needs to be compiled first - see
# cuser_feat.py for suitable commands.
import _cuser_form

env = Environ()

env.io.atom_files_directory = ['../atom_files']
log.verbose()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

class MyGauss(forms.RestraintForm):
    """An implementation of Modeller's harmonic/Gaussian restraint (type 3)
       as a C extension module"""

    _builtin_index = _cuser_form.myform_create()

    def __init__(self, group, feature, mean, stdev):
        forms.RestraintForm.__init__(self, group, feature, 0, (mean, stdev))


mdl = complete_pdb(env, "1fdn")
sel = Selection(mdl)
rsr = mdl.restraints
at = mdl.atoms
rsr.add(MyGauss(group=physical.bond,
                feature=features.Distance(at['CB:1:A'], at['CA:1:A']),
                mean=1.5380, stdev=0.0364))
sel.energy()
