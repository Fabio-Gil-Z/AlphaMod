import _modeller
from modeller.util import modutil
from modeller.util.deprecation import _deprecation_handler


class Feature(object):
    """Base class for all features defined on model atoms"""

    _builtin_index = None
    _user_index = None
    numatoms = 0

    def __init__(self, *atoms):
        if '_builtin_index' not in self.__class__.__dict__ \
           and self._user_index is None:
            self.__class__._user_index = \
                _modeller.mod_user_feature_new(self.eval, self.deriv,
                                               self.is_angle)
        (self.__atoms, self.__mdl) = self.__get_list_atom_indices(atoms)
        self._check_atoms(self.__atoms)

    def get_value(self, deriv=False):
        """Get the current feature value, and also derivatives (with respect
           to atomic coordinates) if ``deriv`` is True."""
        (inds, mdl) = self.get_atom_indices()
        val = self.eval(mdl, inds)
        if deriv:
            return (val, self.deriv(mdl, inds, val))
        else:
            return val

    def is_angle(self):
        """Return True if the feature is an angle (i.e. measures in radians,
           and has a 2pi periodicity)."""
        if self._builtin_index is not None:
            return _modeller.mod_feature_isangle(self._builtin_index)
        else:
            raise NotImplementedError("is_angle not defined for feature "
                                      + str(self))

    def get_atom_indices(self):
        """Returns the atom indices this feature is defined on"""
        return (self.__atoms, self.__mdl)

    def __get_list_atom_indices(self, atoms):
        inds = []
        mdl = None
        for obj in atoms:
            if modutil.non_string_iterable(obj) \
               and not hasattr(obj, 'get_atom_indices'):
                (objinds, objmdl) = self.__get_list_atom_indices(obj)
            else:
                (objinds, objmdl) = obj.get_atom_indices()
            if mdl is None:
                mdl = objmdl
            elif mdl != objmdl:
                raise ValueError("All atoms must be from the same model")
            inds.extend(objinds)
        return (inds, mdl)

    def _check_atoms(self, atoms):
        if len(atoms) != self.numatoms:
            raise ValueError("This feature is defined on %d atoms - got %d" %
                             (self.numatoms, len(atoms)))

    def get_type(cls):
        """Get the numeric type of this feature"""
        if cls._user_index is not None:
            return cls._user_index
        elif cls._builtin_index is not None:
            return cls._builtin_index
        else:
            raise ValueError("Cannot get type - no objects created yet")
    get_type = classmethod(get_type)

    def eval(self, mdl, atom_indices):
        if self._builtin_index is not None:
            return _modeller.mod_feature_eval(mdl.modpt, atom_indices,
                                              self._builtin_index)
        else:
            raise NotImplementedError("Value not defined for feature "
                                      + str(self))

    def deriv(self, mdl, atom_indices, feat):
        if self._builtin_index is not None:
            return _modeller.mod_feature_deriv(mdl.modpt, atom_indices,
                                               self._builtin_index, feat)
        else:
            raise NotImplementedError("Derivatives not defined for feature "
                                      + str(self))

    def indices_to_atoms(self, mdl, atom_indices):
        """Converts Modeller-style atom indices into atom objects"""
        return [mdl.atoms[x-1] for x in atom_indices]


class Distance(Feature):
    """Cartesian distance between two atoms"""
    _builtin_index = 1
    numatoms = 2


class Angle(Feature):
    """Angle (in radians) between three atoms"""
    _builtin_index = 2
    numatoms = 3


class Dihedral(Feature):
    """Dihedral angle (in radians) between four atoms"""
    _builtin_index = 3
    numatoms = 4


class MinimalDistance(Feature):
    """Distance between the nearest pair of atoms"""
    _builtin_index = 6

    def _check_atoms(self, atoms):
        if len(atoms) % 2 != 0 or len(atoms) == 0:
            raise ValueError("This feature is defined on a non-zero even " +
                             "number of atoms - got %d" % len(atoms))


class SolventAccess(Feature):
    """Atomic area exposed to solvent"""
    _builtin_index = 7
    numatoms = 1


class Density(Feature):
    """Atomic density (number of atoms within contact_shell)"""
    _builtin_index = 8
    numatoms = 1


class XCoordinate(Feature):
    """x coordinate of an atom"""
    _builtin_index = 9
    numatoms = 1


class YCoordinate(Feature):
    """y coordinate of an atom"""
    _builtin_index = 10
    numatoms = 1


class ZCoordinate(Feature):
    """z coordinate of an atom"""
    _builtin_index = 11
    numatoms = 1


class DihedralDiff(Feature):
    """Difference (in radians) between two dihedrals"""
    _builtin_index = 12
    numatoms = 8


# Modeller 9 compatibility
class feature(Feature):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(feature)
        Feature.__init__(self, *args, **keys)


class distance(Distance):
    _builtin_index = Distance._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(distance)
        Distance.__init__(self, *args, **keys)


class angle(Angle):
    _builtin_index = Angle._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(angle)
        Angle.__init__(self, *args, **keys)


class dihedral(Dihedral):
    _builtin_index = Dihedral._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(dihedral)
        Dihedral.__init__(self, *args, **keys)


class minimal_distance(MinimalDistance):
    _builtin_index = MinimalDistance._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(minimal_distance)
        MinimalDistance.__init__(self, *args, **keys)


class solvent_access(SolventAccess):
    _builtin_index = SolventAccess._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(solvent_access)
        SolventAccess.__init__(self, *args, **keys)


class density(Density):
    _builtin_index = Density._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(density)
        Density.__init__(self, *args, **keys)


class x_coordinate(XCoordinate):
    _builtin_index = XCoordinate._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(x_coordinate)
        XCoordinate.__init__(self, *args, **keys)


class y_coordinate(YCoordinate):
    _builtin_index = YCoordinate._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(y_coordinate)
        YCoordinate.__init__(self, *args, **keys)


class z_coordinate(ZCoordinate):
    _builtin_index = ZCoordinate._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(z_coordinate)
        ZCoordinate.__init__(self, *args, **keys)


class dihedral_diff(DihedralDiff):
    _builtin_index = DihedralDiff._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(dihedral_diff)
        DihedralDiff.__init__(self, *args, **keys)
