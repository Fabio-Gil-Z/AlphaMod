"""Classes for handling rigid bodies"""
from modeller.util.modobject import ModObject
from modeller.util.deprecation import _deprecation_handler

__docformat__ = "epytext en"


class RigidBody(ModObject):
    """A group of atoms to treat as a rigid body"""

    #: Scaling factor from system state (as used in conjugate gradients and
    #: quasi Newton optimization) to orientation in radians
    scale_factor = 1.0
    _atoms = None

    def __init__(self, *atoms):
        self._atoms = atoms

    def _get_base_atoms(self, mdl):
        return mdl.get_list_atom_indices(self._atoms, None)

    def __get_atoms(self):
        return self._atoms
    atoms = property(__get_atoms, doc="Atoms within the rigid body")


# Modeller 9 compatibility
class rigid_body(RigidBody):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(rigid_body)
        RigidBody.__init__(self, *args, **keys)
