import _modeller
from modeller.util.modobject import ModObject
from modeller.util.logger import log
from modeller.util.deprecation import _deprecation_handler


class IOData(ModObject):
    """Controls reading from/writing to atom files"""

    hetatm = False
    hydrogen = False
    water = False
    convert_modres = True
    hybrid36 = True
    two_char_chain = True
    __atom_files_directory = ['']
    __modpt = None
    __free_func = _modeller.mod_io_data_free

    def __init__(self, copy=None, **kwargs):
        self.__modpt = _modeller.mod_io_data_new()
        if copy:
            for member in copy.__dict__:
                if 'IOData' not in member:
                    self.__dict__[member] = copy.__dict__[member]
            # Copy the list, rather than making another reference
            self.atom_files_directory = copy.atom_files_directory[:]
        for key in kwargs:
            if key in dir(IOData):
                exec("self."+key+"="+repr(kwargs[key]))
            else:
                raise KeyError(str(key))

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_io_data_new()

    def __del__(self):
        if self.__modpt:
            self.__free_func(self.__modpt)

    def __get_atfil(self):
        return self.__atom_files_directory

    def __set_atfil(self, val):
        if isinstance(val, tuple):
            self.__atom_files_directory = list(val)
        elif isinstance(val, list):
            self.__atom_files_directory = val
        elif isinstance(val, str):
            spl = val.split(':')
            log.warning("IOData",
"""Setting io.atom_files_directory to a colon-delimited string is
              deprecated, as it is not robust on Windows systems. Set it to
              a list of directories instead. For example:
              env.io.atom_files_directory = """ + str(spl))  # noqa: E128
            self.__atom_files_directory = spl
        else:
            raise TypeError("atom_files_directory should be a list")

    def __get_modpt(self):
        modpt = self.__modpt
        _modeller.mod_io_data_set(modpt, self.hydrogen, self.hetatm,
                                  self.water, self.convert_modres,
                                  self.hybrid36, self.two_char_chain,
                                  self.__atom_files_directory)
        return modpt

    atom_files_directory = property(__get_atfil, __set_atfil,
                                    doc="List of paths to search "
                                        "for atom files")
    modpt = property(__get_modpt)


# Modeller 9 compatibility
class io_data(IOData):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(io_data)
        IOData.__init__(self, *args, **keys)
