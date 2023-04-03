"""Classes to handle information from CHARMM parameter files"""

import _modeller
from modeller import modfile

__docformat__ = "epytext en"


class Parameters(object):
    """All information from a parameter file. You should never need to
       create a :class:`Parameters` object yourself - one is created for
       you by the :class:`Environ` class, e.g. as ``env.libs.parameters``."""

    def __init__(self, libs):
        self.__libs = libs
        self._modpt = _modeller.mod_libraries_prm_get(libs.modpt)

    def clear(self):
        """Clear all parameter information"""
        _modeller.mod_parameters_clear(self._modpt)

    def append(self, file):
        """Append information from a CHARMM parameter file"""
        fh = modfile._get_filehandle(file, 'r')
        return _modeller.mod_parameters_read(self._modpt, self.__libs.modpt,
                                             fh.file_pointer)

    def read(self, file):
        """Read information from a CHARMM parameter file"""
        self.clear()
        return self.append(file)

    def __get_in_memory(self):
        return _modeller.mod_parameters_in_memory(self._modpt)

    in_memory = property(__get_in_memory,
                         doc="True if information has been read into memory")
