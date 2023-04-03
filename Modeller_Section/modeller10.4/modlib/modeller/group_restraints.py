import _modeller
from modeller import modfile
from modeller.util.modlist import FixList
from modeller.util.modobject import ModObject
from modeller.util.deprecation import _deprecation_handler


class _AtomType(object):
    def __init__(self, residue, atom):
        self.residue = residue
        self.atom = atom

    def __repr__(self):
        return "<%s>" % str(self)

    def __str__(self):
        return "%s:%s" % (self.residue, self.atom)


class _AtomClass(object):
    def __init__(self, atclass, indx):
        self.__atclass = atclass
        self.__indx = indx

    def __get_name(self):
        return _modeller.mod_atom_classes_name_get(self.__atclass._modpt,
                                                   self.__indx)

    def __get_atom_types(self):
        atc = self.__atclass._modpt
        iattmod = _modeller.mod_atom_classes_iattmod_get(atc)
        ngrpatm = _modeller.mod_atom_classes_ngrpatm_get(atc)
        first_entry = _modeller.mod_int1_get(iattmod, self.__indx) - 1
        if self.__indx == ngrpatm - 1:
            last_entry = _modeller.mod_atom_classes_nattmod_get(atc)
        else:
            last_entry = _modeller.mod_int1_get(iattmod, self.__indx + 1) - 1
        types = []
        for i in range(first_entry, last_entry):
            types.append(_AtomType(
                              _modeller.mod_atom_classes_residue_get(atc, i),
                              _modeller.mod_atom_classes_atom_get(atc, i)))
        return tuple(types)

    name = property(__get_name, doc="Class name")
    atom_types = property(__get_atom_types, doc="All atom types in this class")


class _AtomClassList(FixList):
    def __init__(self, gprsr):
        # Keep a reference to make sure GroupRestraints aren't freed
        self.__gprsr = gprsr
        self._modpt = _modeller.mod_group_restraints_atclass_get(gprsr._modpt)

    def __len__(self):
        return _modeller.mod_atom_classes_ngrpatm_get(self._modpt)

    def _getfunc(self, indx):
        return _AtomClass(self, indx)


class GroupRestraints(ModObject):
    """Holds restraints which act on atom classes/groups"""
    _modpt = None
    env = None
    __free_func = _modeller.mod_group_restraints_free
    __class_file = None
    __param_files = None

    def __init__(self, env, classes, parameters=None):
        self.env = env.copy()
        self._modpt = _modeller.mod_group_restraints_new(self)
        self.__read_classes(classes)
        self.__param_files = []
        if parameters:
            self.append(parameters)

    def __del__(self):
        if self._modpt:
            self.__free_func(self._modpt)

    def __setstate__(self, d):
        """Restore internal information from files"""
        self.__dict__.update(d)
        self._modpt = _modeller.mod_group_restraints_new(self)
        self.__read_classes(self.__class_file)
        params = self.__param_files
        self.__param_files = []
        for file in params:
            self.append(file)

    def __read_classes(self, file):
        """Read atom classes from a file"""
        self.__class_file = file
        fh = modfile._get_filehandle(file, 'r')
        return _modeller.mod_atom_classes_read(self._modpt, fh.file_pointer)

    def append(self, file):
        """Read interaction parameters from a file"""
        self.__param_files.append(file)
        fh = modfile._get_filehandle(file, 'r')
        return _modeller.mod_group_restraints_read(self._modpt,
                                                   self.env.libs.modpt,
                                                   fh.file_pointer)

    def __get_atclass(self):
        return _AtomClassList(self)
    atom_classes = property(__get_atclass, doc="Atom classes")


# Modeller 9 compatibility
class group_restraints(GroupRestraints):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(group_restraints)
        GroupRestraints.__init__(self, *args, **keys)
