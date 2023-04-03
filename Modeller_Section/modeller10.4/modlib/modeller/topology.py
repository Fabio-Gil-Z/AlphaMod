"""Classes to handle information from residue topology files"""

import _modeller
from modeller import modfile


class Topology(object):
    """All information from a residue topology file. You should never need to
       create a :class:`Topology` object yourself - one is created for you
       by the :class:`Environ` class, e.g. as ``env.libs.topology``."""

    def __init__(self, libs):
        self.__libs = libs
        self._modpt = _modeller.mod_libraries_tpl_get(libs.modpt)

    def clear(self):
        """Remove all topology information"""
        _modeller.mod_topology_clear(self._modpt, self.__libs.modpt)

    def read(self, file):
        """Reads a residue topology file"""
        self.clear()
        return self.append(file)

    def append(self, file):
        """Appends information from a residue topology file"""
        fh = modfile._get_filehandle(file, 'r')
        return _modeller.mod_topology_read(self._modpt, self.__libs.modpt,
                                           fh.file_pointer)

    def make(self, submodel):
        """Reduces the most detailed topology to a sub-topology model.

           :param int submodel: an integer from 1 to 10, to specify the
                  sub-topology model as defined in ``models.lib``.
        """
        self.submodel = submodel
        return _modeller.mod_topology_model_make(self._modpt,
                                                 self.__libs.modpt)

    def write(self, file):
        """Writes topology library to a file"""
        fh = modfile._get_filehandle(file, 'w')
        return _modeller.mod_topology_model_write(self._modpt,
                                                  self.__libs.modpt,
                                                  fh.file_pointer)

    def __get_submodel(self):
        return _modeller.mod_topology_submodel_get(self._modpt)

    def __set_submodel(self, val):
        _modeller.mod_topology_submodel_set(self._modpt, val)

    def __get_in_memory(self):
        return _modeller.mod_topology_in_memory(self._modpt)

    submodel = property(__get_submodel, __set_submodel,
                        doc="Topology sub-model, an integer from 1 to 10, "
                            "as defined in ``models.lib``")
    in_memory = property(__get_in_memory,
                         doc="True if information has been read into memory")
