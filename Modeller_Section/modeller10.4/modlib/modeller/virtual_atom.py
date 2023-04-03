from modeller.pseudo_atom import PseudoAtom
from modeller.util.deprecation import _deprecation_handler


class CH1(PseudoAtom):
    """Virtual aliphatic proton on a tetrahedral carbon (->CH),
       defined by the central C and the three other substituents"""
    _builtin_index = 2
    numatoms = 4


class CH1A(PseudoAtom):
    """Virtual aromatic proton on a trigonal carbon (=CH),
       defined by the central C and the two C atoms bonded to the central C"""
    _builtin_index = 3
    numatoms = 3


class CH2(PseudoAtom):
    """Virtual aliphatic proton on a tetrahedral carbon (>CH2)
       assigned stereospecifically; defined by the central
       tetrahedral atom and the other two substituents on it"""
    _builtin_index = 5
    numatoms = 3


# Modeller 9 compatibility
class ch1(CH1):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(ch1)
        CH1.__init__(self, *args, **keys)


class ch1a(CH1A):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(ch1a)
        CH1A.__init__(self, *args, **keys)


class ch2(CH2):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(ch2)
        CH2.__init__(self, *args, **keys)
