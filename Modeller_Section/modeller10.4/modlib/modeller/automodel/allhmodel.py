"""Classes to build all-atom model(s) using template information"""

from modeller.automodel.automodel import AutoModel
from modeller.util.deprecation import _deprecation_handler

__docformat__ = "epytext en"


class AllHModel(AutoModel):
    """Automatically build all-atom model(s) using template information"""

    toplib = '${LIB}/top_allh.lib'

    def __init__(self, env, alnfile, knowns, sequence, deviation=None,
                 library_schedule=None, csrfile=None, inifile=None,
                 assess_methods=None):
        AutoModel.__init__(self, env, alnfile, knowns, sequence, deviation,
                           library_schedule, csrfile, inifile, assess_methods)
        # Modeling won't work unless we read/write hydrogen atoms!
        self.env.io.hydrogen = True


# Modeller 9 compatibility
class allhmodel(AllHModel):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(allhmodel)
        AllHModel.__init__(self, *args, **keys)
