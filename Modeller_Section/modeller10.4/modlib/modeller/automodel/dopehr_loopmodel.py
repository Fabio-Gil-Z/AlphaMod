"""Classes for loop modeling using the DOPE-HR potential"""

from modeller.automodel.dope_loopmodel import DOPELoopModel
from modeller import GroupRestraints
from modeller.util.deprecation import _deprecation_handler

__docformat__ = "epytext en"


class DOPEHRLoopModel(DOPELoopModel):
    """Loop modeling using the DOPE-HR potential"""

    def read_potential(self):
        return GroupRestraints(self.env, classes='$(LIB)/atmcls-mf.lib',
                               parameters='$(LIB)/dist-mfhr.lib')


# Modeller 9 compatibility
class dopehr_loopmodel(DOPEHRLoopModel):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(dopehr_loopmodel)
        DOPEHRLoopModel.__init__(self, *args, **keys)
