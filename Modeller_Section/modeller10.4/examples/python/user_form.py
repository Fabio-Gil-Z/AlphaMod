from modeller import *
from modeller.scripts import complete_pdb

env = Environ()

env.io.atom_files_directory = ['../atom_files']
log.verbose()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

class MyGauss(forms.RestraintForm):
    """An implementation of Modeller's harmonic/Gaussian restraint (type 3)
       in pure Python"""

    rt = 0.5900991    # RT at 297.15K, in kcal/mol

    def __init__(self, group, feature, mean, stdev):
        forms.RestraintForm.__init__(self, group, feature, 0, (mean, stdev))

    def eval(self, feats, iftyp, modal, param, deriv):
        (mean, stdev) = param
        delt = self.deltaf(feats[0], mean, iftyp[0])
        val = self.rt * 0.5 * delt**2  / stdev**2
        if deriv:
            fderv = self.rt * delt / stdev**2
            return val, [fderv]
        else:
            return val

    def vmin(self, feats, iftyp, modal, param):
        (mean, stdev) = param
        return self.deltaf(feats[0], mean, iftyp[0])

    def rvmin(self, feats, iftyp, modal, param):
        (mean, stdev) = param
        return self.deltaf(feats[0], mean, iftyp[0]) / stdev

    def min_mean(self, feats, iftyp, modal, param):
        (mean, stdev) = param
        return [mean]

    def get_range(self, iftyp, modal, param, spline_range):
        (mean, stdev) = param
        return (mean - stdev * spline_range, mean + stdev * spline_range)

    # There is only one minimum, so the 'heavy' mean is the same as the 'min'
    vheavy = vmin
    rvheavy = rvmin
    heavy_mean = min_mean

mdl = complete_pdb(env, "1fdn")
sel = Selection(mdl)
rsr = mdl.restraints
at = mdl.atoms
rsr.add(MyGauss(group=physical.bond,
                feature=features.Distance(at['CB:1:A'], at['CA:1:A']),
                mean=1.5380, stdev=0.0364))
sel.energy()

# Restraints using user-defined forms can be converted to splines for speed.
# Only one-dimensional forms that define the get_range() method can be splined.
rsr.spline(MyGauss, features.Distance, physical.bond, spline_dx=0.05)
sel.energy()
