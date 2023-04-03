/* Convert Python strings to/from C fixed-size buffers */


#ifdef SWIGPYTHON
%typemap(in) (const void *buffer, size_t bufsize) (Py_ssize_t temp) {
%#if PY_VERSION_HEX < 0x03000000
  if (PyString_AsStringAndSize($input, (char **)&$1, &temp) == -1) {
    SWIG_fail;
  } else {
    $2 = (size_t)temp;
  }
%#else
  if (PyBytes_AsStringAndSize($input, (char **)&$1, &temp) == -1) {
    SWIG_fail;
  } else {
    $2 = (size_t)temp;
  }
%#endif
}

%typemap(in,fragment=SWIG_AsVal_frag(size_t)) (void *buffer, size_t bufsize) (int res) {
  res = SWIG_AsVal(size_t)($input, &$2);
  if (!SWIG_IsOK(res)) {
    %argument_fail(res, "(size_t bufsize)", $symname, $argnum);
  }
  $1 = g_malloc($2);
}

%typemap(freearg,match="in") (void *buffer, size_t bufsize) {
  g_free($1);
}

%typemap(argout) (void *buffer, size_t bufsize) {
%#if PY_VERSION_HEX < 0x03000000
  %append_output(PyString_FromStringAndSize((char *)$1, (Py_ssize_t)$2));
%#else
  %append_output(PyBytes_FromStringAndSize((char *)$1, (Py_ssize_t)$2));
%#endif
}
#endif
