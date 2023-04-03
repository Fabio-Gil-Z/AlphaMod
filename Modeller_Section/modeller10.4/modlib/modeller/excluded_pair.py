from modeller.util.deprecation import _deprecation_handler


class ExcludedPair(object):
    """A pair of atoms to exclude from the non-bonded list"""

    def __init__(self, *atoms):
        self._atoms = atoms

    def _get_base_atoms(self, mdl):
        return mdl.get_list_atom_indices(self._atoms, 2)

    def __get_atoms(self):
        return self._atoms
    atoms = property(__get_atoms, doc="Excluded atoms")


# Modeller 9 compatibility
class excluded_pair(ExcludedPair):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(excluded_pair)
        ExcludedPair.__init__(self, *args, **keys)
