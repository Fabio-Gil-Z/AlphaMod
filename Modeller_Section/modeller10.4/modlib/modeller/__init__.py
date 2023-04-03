"""MODELLER, a package for protein structure modeling

   See U{https://salilab.org/modeller/} for further details.

   @author: Andrej Sali
   @copyright: 1989-2022 Andrej Sali
"""

__all__ = ['EnergyData', 'IOData', 'GroupRestraints', 'Symmetry', 'RigidBody',
           'ExcludedPair', 'ModellerError', 'FileFormatError',
           'StatisticsError', 'SequenceMismatchError', 'Model', 'Alignment',
           'Environ', 'SequenceDB', 'Profile', 'SAXSData', 'Density', 'PSSMDB',
           'Selection', 'log', 'info', 'pseudo_atom', 'virtual_atom',
           'modfile', 'features', 'forms', 'secondary_structure', 'terms',
           'physical']

__docformat__ = "epytext en"

# Version check
import sys
if sys.version_info[0] < 2 \
   or (sys.version_info[0] == 2 and sys.version_info[1] < 3):
    raise ImportError("This module requires Python 2.3 or later")

try:
    from modeller import config
except ImportError:
    config = None

def __get_python_api_ver():
    """Get the Python version at which the C API last changed."""
    ver = sys.version_info[0:2]
    if ver[0] == 3:
        if ver[1] < 2:
            return (3,0)
        elif ver[1] == 2:
            return (3,2)
        else:
            return (3,3)
    elif ver[0] == 2 and ver[1] >= 5:
        return (2,5)
    return None

def __is_64_bit_windows():
    """Return True iff we are running on 64-bit Windows."""
    # This only works on Python 2.5 or later
    if hasattr(sys, 'maxsize'):
        return sys.maxsize > 2**32
    # This works on older Pythons, but not in Python 3
    else:
        return type(sys.dllhandle) == long

# Special processing on Windows to find _modeller.pyd and Modeller DLLs:
if hasattr(config, 'install_dir') and hasattr(sys, 'dllhandle'):
    if __is_64_bit_windows():
        exetype = 'x86_64-w64'
    else:
        exetype = 'i386-w32'
    dpath = config.install_dir + '\\lib\\%s\\python%d.%d' \
                                 % ((exetype,) + sys.version_info[:2])
    if dpath not in sys.path:
        # Insert *after* first entry, so as not to break parallel module's
        # propagation of the first entry (that containing the running script's
        # directory) to workers
        sys.path.insert(1, dpath)
    try:
        import os
        dpath = config.install_dir + '\\lib\\%s' % exetype
        if dpath not in os.environ['PATH']:
            os.environ['PATH'] = dpath + ';' + os.environ['PATH']
        # Python 3.8 or later don't look in PATH for DLLs
        if hasattr(os, 'add_dll_directory'):
            __dll_directory = os.add_dll_directory(dpath)
        del os
    except ImportError:
        pass
    del dpath
# Add Python version-specific directory to search path:
elif hasattr(config, 'install_dir'):
    api_ver = __get_python_api_ver()
    if api_ver is not None:
        try:
            import os.path, re, sys
            srch = re.compile("%s/*lib/[^/]+/?" % config.install_dir)
            for (n, pathcomp) in enumerate(sys.path):
                if srch.match(pathcomp):
                    modpath = os.path.join(pathcomp, 'python%d.%d' % api_ver)
                    if modpath not in sys.path and os.path.exists(modpath):
                        sys.path.insert(n, modpath)
                    break
            del re, n, pathcomp, os, srch
        except ImportError:
            pass

# Set Modeller install location and license
import _modeller
if hasattr(config, 'license'):
    _modeller.mod_license_key_set(config.license)
if hasattr(config, 'install_dir'):
    _modeller.mod_install_dir_set(config.install_dir)

_modeller.mod_start()
__version__ = _modeller.mod_short_version_get()

from modeller.energy_data import EnergyData
from modeller.io_data import IOData
from modeller.environ import Environ
from modeller.group_restraints import GroupRestraints
from modeller.error import *
from modeller.model import Model
from modeller.alignment import Alignment
from modeller.sequence_db import SequenceDB
from modeller.profile import Profile
from modeller.saxsdata import SAXSData
from modeller.density import Density
from modeller.pssmdb import PSSMDB
from modeller.excluded_pair import ExcludedPair
from modeller.rigid_body import RigidBody
from modeller.symmetry import Symmetry
from modeller.selection import Selection
from modeller.util.logger import log
from modeller.information import info
from modeller import pseudo_atom
from modeller import virtual_atom
from modeller import modfile
from modeller import features
from modeller import forms
from modeller import secondary_structure
from modeller import terms
from modeller import physical

# Modeller 9 compatibility
from modeller.io_data import io_data
from modeller.energy_data import energy_data
from modeller.symmetry import symmetry
from modeller.rigid_body import rigid_body
from modeller.excluded_pair import excluded_pair
from modeller.group_restraints import group_restraints
from modeller.density import density
from modeller.saxsdata import saxsdata
from modeller.pssmdb import pssmdb
from modeller.profile import profile
from modeller.sequence_db import sequence_db
from modeller.selection import selection
from modeller.alignment import alignment
from modeller.model import model
from modeller.environ import environ
__all__.extend(['io_data', 'energy_data', 'symmetry', 'rigid_body',
                'excluded_pair', 'group_restraints', 'density', 'saxsdata',
                'pssmdb', 'profile', 'sequence_db', 'selection', 'alignment',
                'model', 'environ'])

# Load in readline, if available, to make interactive use easier
try:
    if len(sys.argv) > 0 and sys.argv[0] == '-' and sys.stdin.isatty():
        import readline
except (ImportError, AttributeError):
    pass

# Set job name
if len(sys.argv) > 0 and sys.argv[0] != '-':
    nam = sys.argv[0]
    if nam.endswith('.py'):
        nam = nam[:-3]
    info.jobname = nam
    del nam

del sys, _modeller, config
