/* Conversions between Python strings and C GString objects */

%typemap(in) GString * (PyObject *temp) {
  $1 = NULL;
%#if PY_VERSION_HEX < 0x03000000
  if (PyString_Check($input)) {
    temp = $input;
  } else {
    PyErr_Format(PyExc_TypeError, "Expected a string for argument number %d",
                 $argnum);
    SWIG_fail;
  }
%#else
  if (PyUnicode_Check($input)) {
    temp = PyUnicode_AsUTF8String($input);
    if (!temp) {
      SWIG_fail;
    }
  } else {
    PyErr_Format(PyExc_TypeError, "Expected a string for argument number %d",
                 $argnum);
    SWIG_fail;
  }
%#endif
}

%typemap(memberin) GString * {
%#if PY_VERSION_HEX < 0x03000000
  g_string_assign($1, PyString_AsString(temp2));
%#else
  g_string_assign($1, PyBytes_AsString(temp2));
  Py_DECREF(temp2);
%#endif
}

%typemap(out) GString * {
%#if PY_VERSION_HEX < 0x03000000
  $result = PyString_FromString($1->str);
%#else
  $result = PyUnicode_FromString($1->str);
%#endif
}
