import _modeller
from modeller.util.deprecation import _deprecation_handler


class Symmetry(object):
    """Symmetry restraint"""

    def __init__(self, set1, set2, weight):
        self.__segments = []
        self.append(set1, set2, weight)

    def append(self, set1, set2, weight):
        """Add a pair of atom sets to the restraint"""
        self.__segments.append((set1, set2, weight))

    def _add_segments(self, mdl):
        """Internal function to add the segments to the given model's
           restraints list"""
        addtoseg = False
        for (set1, set2, weight) in self.__segments:
            inds1 = mdl.get_list_atom_indices(set1, None)
            inds2 = mdl.get_list_atom_indices(set2, None)
            _modeller.mod_symmetry_define(mdl.modpt, inds1, inds2,
                                          (True, addtoseg), weight)
            addtoseg = True


# Modeller 9 compatibility
class symmetry(Symmetry):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(symmetry)
        Symmetry.__init__(self, *args, **keys)
