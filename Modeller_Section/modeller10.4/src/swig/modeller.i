%module _modeller

/* '#' formats in parsing or building Python values
   (e.g. in PyObject_CallFunction) use Py_ssize_t, not int, for lengths
   in Python >= 2.5 */
%begin %{
#define PY_SSIZE_T_CLEAN
%}

%{
#include "modeller.h"
#if PY_VERSION_HEX < 0x03000000
# include "cStringIO.h"
#endif
%}

%{
#ifdef SWIGPYTHON
#if PY_VERSION_HEX >= 0x03000000
/* The io.TextIOBase class */
static PyObject *text_io_base = NULL;
#endif
#endif
%}

%include "custom-typemaps.i"

typedef int gboolean;

%ignore mod_energy_term_add_profile;

/* Wrap some Modeller functions and structures directly from the header files */
%include "../include/modeller.h"
%include "../include/mod_core.h"
%include "../include/mod_selection.h"
%include "../include/mod_libs.h"
%include "../include/mod_build.h"
%include "../include/mod_info.h"
%include "../include/mod_system.h"
%include "../include/mod_ga341.h"
%include "../include/mod_user_features.h"
%include "../include/mod_user_forms.h"
%include "../include/mod_user_terms.h"
%include "../include/mod_user_simloc.h"
%include "../include/mod_user_distsimloc.h"
%include "../include/mod_optimizer_actions.h"
%include "../include/mod_gbsa.h"
%include "../include/mod_soap_od.h"
%include "../include/mod_soap_access.h"
%include "../include/mod_soap_pair.h"
%include "../include/mod_libraries.h"

/* Provide custom wrapping of some Modeller structures and functions */
%include "mod_log.i"
%include "mod_file.i"
%include "optimizer.i"
%include "coordinates.i"
%include "model.i"
%include "alignment.i"
%include "profile.i"
%include "topology.i"
%include "density.i"
%include "sequence.i"
%include "alnsequence.i"
%include "sequence_db.i"
%include "sequence_db_window.i"
%include "saxsdata.i"
%include "pseudo_atoms.i"
%include "symmetry.i"
%include "restraints.i"
%include "fortran-pointers.i"
%include "model_topology.i"
%include "structure.i"
%include "libraries.i"
%include "parameters.i"
%include "schedule.i"
%include "energy_data.i"
%include "io_data.i"
%include "group_restraints.i"

%init {
#ifdef SWIGPYTHON
%#if PY_VERSION_HEX < 0x03000000
  /* Initialize cStringIO API */
  PycString_IMPORT;
  /* If not available, ignore the ImportError - we'll handle it later */
  if (PyErr_Occurred()) {
    PyErr_Clear();
  }
%#else
  /* Get io.TextIOBase class */
  {
    PyObject *io = PyImport_ImportModule("io");
    if (io) {
      text_io_base = PyObject_GetAttrString(io, "TextIOBase");
      /* todo: decref text_io_base when Modeller is unloaded */
      Py_DECREF(io);
    }
  }
%#endif

  /* Create base error class */
  moderror = PyErr_NewException("_modeller.ModellerError", NULL, NULL);
  Py_INCREF(moderror);
  PyModule_AddObject(m, "ModellerError", moderror);

  /* Create error subclasses */
  file_format_error = PyErr_NewException("_modeller.FileFormatError", moderror,
                                         NULL);
  Py_INCREF(file_format_error);
  PyModule_AddObject(m, "FileFormatError", file_format_error);

  statistics_error = PyErr_NewException("_modeller.StatisticsError", moderror,
                                        NULL);
  Py_INCREF(statistics_error);
  PyModule_AddObject(m, "StatisticsError", statistics_error);

  sequence_mismatch_error
      = PyErr_NewException("_modeller.SequenceMismatchError", moderror, NULL);
  Py_INCREF(sequence_mismatch_error);
  PyModule_AddObject(m, "SequenceMismatchError", sequence_mismatch_error);
#endif
}
