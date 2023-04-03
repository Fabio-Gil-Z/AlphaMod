"""Classes for obtaining per-residue energy profiles."""

from modeller.util.modlist import FixList


def __smooth_one(vals, n, window):
    """Return the weighted average of ``window`` values around ``vals[n]``"""
    i1 = max(0, n - window)
    i2 = min(len(vals) - 1, n + window)
    wtsum = av = 0.
    for j in range(i1, i2 + 1):
        weight = 0.1 * (window - abs(j - n) + 1)
        av += weight * vals[j]
        wtsum += weight
    return av / max(wtsum, 0.1)


def _smooth(vals, window):
    """Return a smoothed copy of the array ``vals``"""
    return [__smooth_one(vals, n, window) for n in range(len(vals))]


class EnergyProfileResidue(object):
    """A single residue in the energy profile."""
    def __init__(self, energy, num_restraints):
        self.energy = energy
        self.num_restraints = num_restraints

    def __repr__(self):
        return "<%6.3f / %d>" % (self.energy, self.num_restraints)

    def __get_normalized_energy(self):
        if self.num_restraints == 0:
            return 0.
        else:
            return self.energy / self.num_restraints
    normalized_energy = property(__get_normalized_energy,
                                 doc="'Normalized' energy, where the energy"
                                     " is divided by the number of"
                                     " restraints acting on this residue")


class EnergyProfile(FixList):
    """A per-residue energy profile"""

    def __init__(self, profile, nprofile, min_rms, heavy_rms):
        self.__profile = [EnergyProfileResidue(energy, num_restraints)
                          for energy, num_restraints in zip(profile, nprofile)]
        self.min_rms = min_rms
        self.heavy_rms = heavy_rms

    def _getfunc(self, indx):
        return self.__profile[indx]

    def __len__(self):
        return len(self.__profile)

    def get_normalized(self):
        """Return a new 'normalized' energy profile, in which each residue's
           energy is divided by the number of restraints acting on that
           residue."""
        return EnergyProfile([r.normalized_energy for r in self.__profile],
                             [1]*len(self), self.min_rms, self.heavy_rms)

    def get_smoothed(self, window=1):
        """Return a new energy profile, smoothed by window averaging."""
        return EnergyProfile(_smooth([r.energy for r in self.__profile],
                                     window),
                             [r.num_restraints for r in self.__profile],
                             self.min_rms, self.heavy_rms)

    def write_to_file(self, filename):
        """Write the profile to a file"""
        if hasattr(filename, 'write'):
            fh = filename
        else:
            fh = open(filename, "w")
        for (n, res) in enumerate(self.__profile):
            fh.write("%10d %12.4f\n" % (n+1, res.energy))
        if not hasattr(filename, 'write'):
            fh.close()
