"""All physical restraint types, used for creating restraints, and in the
   schedule."""

from modeller.util.deprecation import _deprecation_handler


class PhysicalType(object):
    """Physical restraint type"""
    _num = 0
    _types = []
    _shortname = ""

    def __init__(self, doc, shortname):
        PhysicalType._num += 1
        self._num = PhysicalType._num
        self.__doc__ = doc
        self._shortname = shortname
        PhysicalType._types.append(self)

    def get_type(self):
        """Get the numeric type of this physical restraint type"""
        return self._num

    def __str__(self):
        return '<physical restraint type: %s>' % self.__doc__

    def __repr__(self):
        return '<%s>' % self.__doc__
    modpt = property(get_type)


class Values(object):
    """Scaling or cutoff values for all physical restraint types"""

    def __init__(self, default=1.0, **keys):
        self._default = default
        self._dict = {}
        for (term, val) in keys.items():
            term = eval("%s" % term)
            self[term] = val

    def set_from_list(self, lval):
        """Set the state from a Modeller internal list of values"""
        for (num, val) in enumerate(lval):
            typ = PhysicalType._types[num]
            self[typ] = val

    def __len__(self):
        return len(PhysicalType._types)

    def __repr__(self):
        if self._default is not None:
            terms = ["default=%f" % self._default]
        else:
            terms = []
        for term in PhysicalType._types:
            num = term.get_type()
            if num in self._dict:
                terms.append("%s=%f" % (term._shortname, self._dict[num]))
        return "physical.Values(" + ", ".join(terms) + ")"

    def keys(self):
        keys = []
        for term in PhysicalType._types:
            if term.get_type() in self._dict:
                keys.append(term)
        return keys

    def __mul__(self, other):
        obj = Values()
        if not isinstance(other, Values):
            raise TypeError("Must use physical.Values objects here")
        if self._default is None or other._default is None:
            obj._default = None
        else:
            obj._default = self._default * other._default
        for x in self._dict:
            if x in other._dict:
                oval = other._dict[x]
            else:
                if other._default is None:
                    raise ValueError("Default is None, and no value given "
                                     "for %s" % PhysicalType._types[x-1])
                oval = other._default
            obj._dict[x] = self._dict[x] * oval
        for x in other._dict:
            if x not in self._dict:
                if self._default is None:
                    raise ValueError("Default is None, and no value given "
                                     "for %s" % PhysicalType._types[x-1])
                obj._dict[x] = self._default * other._dict[x]
        return obj

    def _typecheck(self, item):
        if not isinstance(item, PhysicalType):
            raise TypeError("keys must be PhysicalType objects or 'default'")

    def __contains__(self, item):
        if isinstance(item, str) and item == 'default':
            return True
        self._typecheck(item)
        return item in PhysicalType._types

    def __getitem__(self, key):
        # Allow read-only access as a list, for the convenience of old code
        # which expects a vector of numbers
        if isinstance(key, int):
            return self[PhysicalType._types[key]]
        if isinstance(key, str) and key == 'default':
            return self._default
        self._typecheck(key)
        try:
            return self._dict[key.get_type()]
        except KeyError:
            return self._default

    def __setitem__(self, key, value):
        if isinstance(key, str) and key == 'default':
            self._default = value
        else:
            self._typecheck(key)
            self._dict[key.get_type()] = value


def from_list(lval):
    """Initialize an object from a Modeller internal list of values"""
    obj = Values(default=None)
    obj.set_from_list(lval)
    return obj


bond = PhysicalType("Bond length potential", "bond")
angle = PhysicalType("Bond angle potential", "angle")
dihedral = PhysicalType("Stereochemical cosine torsion potential", "dihedral")
improper = PhysicalType("Stereochemical improper torsion potential",
                        "improper")
soft_sphere = PhysicalType("Soft-sphere overlap restraints", "soft_sphere")
lennard_jones = PhysicalType("Lennard-Jones 6-12 potential", "lennard_jones")
coulomb = PhysicalType("Coulomb point-point electrostatic potential",
                       "coulomb")
h_bond = PhysicalType("H-bonding potential", "h_bond")
ca_distance = PhysicalType("Distance restraints 1 (CA-CA)", "ca_distance")
n_o_distance = PhysicalType("Distance restraints 2 (N-O)", "n_o_distance")
phi_dihedral = PhysicalType("Mainchain Phi dihedral restraints",
                            "phi_dihedral")
psi_dihedral = PhysicalType("Mainchain Psi dihedral restraints",
                            "psi_dihedral")
omega_dihedral = PhysicalType("Mainchain Omega dihedral restraints",
                              "omega_dihedral")
chi1_dihedral = PhysicalType("Sidechain Chi_1 dihedral restraints",
                             "chi1_dihedral")
chi2_dihedral = PhysicalType("Sidechain Chi_2 dihedral restraints",
                             "chi2_dihedral")
chi3_dihedral = PhysicalType("Sidechain Chi_3 dihedral restraints",
                             "chi3_dihedral")
chi4_dihedral = PhysicalType("Sidechain Chi_4 dihedral restraints",
                             "chi4_dihedral")
disulfide_distance = PhysicalType("Disulfide distance restraints",
                                  "disulfide_distance")
disulfide_angle = PhysicalType("Disulfide angle restraints", "disulfide_angle")
disulfide_dihedral = PhysicalType("Disulfide dihedral angle restraints",
                                  "disulfide_dihedral")
lower_distance = PhysicalType("Lower bound distance restraints",
                              "lower_distance")
upper_distance = PhysicalType("Upper bound distance restraints",
                              "upper_distance")
sd_mn_distance = PhysicalType("Distance restraints 3 (SDCH-MNCH)",
                              "sd_mn_distance")
chi5_dihedral = PhysicalType("Sidechain Chi_5 dihedral restraints",
                             "chi5_dihedral")
phi_psi_dihedral = PhysicalType("Phi/Psi pair of dihedral restraints",
                                "phi_psi_dihedral")
sd_sd_distance = PhysicalType("Distance restraints 4 (SDCH-SDCH)",
                              "sd_sd_distance")
xy_distance = PhysicalType("Distance restraints 5 (X-Y)", "xy_distance")
nmr_distance = PhysicalType("NMR distance restraints 6 (X-Y)", "nmr_distance")
nmr_distance2 = PhysicalType("NMR distance restraints 7 (X-Y)",
                             "nmr_distance2")
min_distance = PhysicalType("Minimal distance restraints", "min_distance")
nonbond_spline = PhysicalType("Non-bonded spline restraints", "nonbond_spline")
accessibility = PhysicalType("Atomic accessibility restraints",
                             "accessibility")
density = PhysicalType("Atomic density restraints", "density")
absposition = PhysicalType("Absolute position restraints", "absposition")
dihedral_diff = PhysicalType("Dihedral angle difference restraints",
                             "dihedral_diff")
gbsa = PhysicalType("GBSA implicit solvent potential", "gbsa")
em_density = PhysicalType("EM density fitting potential", "em_density")
saxs = PhysicalType("SAXS restraints", "saxs")
symmetry = PhysicalType("Symmetry restraints", "symmetry")


# Modeller 9 compatibility
class physical_type(PhysicalType):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(physical_type)
        PhysicalType.__init__(self, *args, **keys)


class values(Values):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(values)
        Values.__init__(self, *args, **keys)
