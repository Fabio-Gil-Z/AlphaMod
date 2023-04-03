"""Classes to cover most comparative modeling tasks.
     - Use the :class:`AutoModel` class to build one or more comparative models.
     - Use the :class:`AllHModel` class instead if you want to build
       all-atom models.
     - Use the :class:`LoopModel`, :class:`DOPELoopModel`, or
       :class:`DOPEHRLoopModel` classes if you additionally want to
       refine loop regions.
"""

from modeller.automodel.automodel import AutoModel
from modeller.automodel.allhmodel import AllHModel
from modeller.automodel.loopmodel import LoopModel
from modeller.automodel.dope_loopmodel import DOPELoopModel
from modeller.automodel.dopehr_loopmodel import DOPEHRLoopModel
from modeller.automodel import refine, generate, randomize
from modeller.automodel import assess, autosched

# Modeller 9 compatibility
from modeller.automodel.automodel import automodel
from modeller.automodel.loopmodel import loopmodel
from modeller.automodel.allhmodel import allhmodel
from modeller.automodel.dope_loopmodel import dope_loopmodel
from modeller.automodel.dopehr_loopmodel import dopehr_loopmodel
