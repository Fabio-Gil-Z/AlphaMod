/* Special handling for exceptions raised by mod_file_open */
%exception mod_file_open {
  $action
  if (!result && check_for_error()) {
    SWIG_fail;
  }
}

/* Don't wrap the mod_file structure; should be treated as opaque */
%ignore mod_file;

/* Provide our own implementation for mod_file_open_stream */
%ignore mod_file_open_stream;

%{
#if PY_VERSION_HEX < 0x03000000
/* Adapter to use cStringIO read with mod_file objects */
static size_t cstr_readfunc(void *data, void *buffer, size_t bufsize, int *ierr)
{
  char *buf;
  int len;
  PyObject *obj = (PyObject *)data;
  len = PycStringIO->cread(obj, &buf, (Py_ssize_t)bufsize);
  if (PyErr_Occurred()) {
    *ierr = 1;
    return 0;
  } else {
    *ierr = 0;
    memcpy(buffer, buf, MIN(len, bufsize));
  }
  return (size_t)len;
}

/* Adapter to use cStringIO readline with mod_file objects */
static void cstr_readlinefunc(void *data, GString *str, int *ierr)
{
  char *buf;
  int len;
  PyObject *obj = (PyObject *)data;
  len = PycStringIO->creadline(obj, &buf);
  if (PyErr_Occurred()) {
    *ierr = 1;
  } else {
    *ierr = 0;
    g_string_truncate(str, 0);
    g_string_append_len(str, buf, len);
  }
}

/* Adapter to use cStringIO write with mod_file objects */
static void cstr_writefunc(void *data, const char *buffer, size_t bufsize,
                           int *ierr)
{
  PyObject *obj = (PyObject *)data;
  PycStringIO->cwrite(obj, (char *)buffer, (Py_ssize_t)bufsize);
  *ierr = (PyErr_Occurred() != NULL);
}

/* Release the reference to the underlying cStringIO object */
static void cstr_freefunc(void *data)
{
  Py_DECREF((PyObject *)data);
}
#endif

struct FileMethods {
  PyObject *read_method, *readline_method, *write_method;
};

static size_t pyfile_readfunc(void *data, void *buffer, size_t bufsize,
                              int *ierr)
{
  struct FileMethods *methods = (struct FileMethods *)data;
  static char fmt[] = "(i)";

  PyObject *result = PyObject_CallFunction(methods->read_method, fmt,
                                           (int)bufsize);
  if (!result) {
    *ierr = 1;
    return 0;
#if PY_VERSION_HEX < 0x03000000
  } else if (PyString_Check(result)) {
    Py_ssize_t len = PyString_GET_SIZE(result);
    memcpy(buffer, PyString_AS_STRING(result), MIN(bufsize, len));
#else
  } else if (PyBytes_Check(result)) {
    Py_ssize_t len = PyBytes_Size(result);
    memcpy(buffer, PyBytes_AS_STRING(result), MIN(bufsize, len));
#endif
    Py_DECREF(result);
    *ierr = 0;
    return (size_t)len;
  } else {
    Py_DECREF(result);
    PyErr_SetString(PyExc_TypeError, "Python file-like object read method "
#if PY_VERSION_HEX < 0x03000000
                    "should return a 'str' object");
#else
                    "should return a 'bytes' object");
#endif
    *ierr = 1;
    return 0;
  }
}

static void pyfile_readlinefunc(void *data, GString *str, int *ierr)
{
  struct FileMethods *methods = (struct FileMethods *)data;
  static char fmt[] = "()";

  PyObject *result = PyObject_CallFunction(methods->readline_method, fmt);
  if (!result) {
    *ierr = 1;
#if PY_VERSION_HEX < 0x03000000
  } else if (PyString_Check(result)) {
    g_string_assign(str, PyString_AS_STRING(result));
#else
  } else if (PyBytes_Check(result)) {
    g_string_assign(str, PyBytes_AsString(result));
#endif
    *ierr = 0;
    Py_DECREF(result);
  } else {
    Py_DECREF(result);
    PyErr_SetString(PyExc_TypeError, "Python file-like object readline method "
#if PY_VERSION_HEX < 0x03000000
                    "should return a 'str' object");
#else
                    "should return a 'bytes' object");
#endif
    *ierr = 1;
  }
}

static void pyfile_writefunc(void *data, const char *buffer, size_t bufsize,
                             int *ierr)
{
  /* If an exception already occurred, exit immediately */
  if (PyErr_Occurred()) {
    *ierr = 1;
  } else {
    struct FileMethods *methods = (struct FileMethods *)data;
#if PY_VERSION_HEX < 0x03000000
    static char fmt[] = "(s#)";
#else
    static char fmt[] = "(y#)";
#endif

    PyObject *result = PyObject_CallFunction(methods->write_method, fmt,
                                             buffer,
#if PY_VERSION_HEX >= 0x02050000 && defined(PY_SSIZE_T_CLEAN)
                                             (Py_ssize_t)bufsize);
#else
                                             (int)bufsize);
#endif
    if (!result) {
      *ierr = 1;
    } else {
      *ierr = 0;
      Py_DECREF(result);
    }
  }
}

#if PY_VERSION_HEX >= 0x03000000
static size_t pyfile_strreadfunc(void *data, void *buffer, size_t bufsize,
                                 int *ierr)
{
  struct FileMethods *methods = (struct FileMethods *)data;
  static char fmt[] = "(i)";

  PyObject *result = PyObject_CallFunction(methods->read_method, fmt,
                                           (int)bufsize);
  if (!result) {
    *ierr = 1;
    return 0;
  } else if (PyUnicode_Check(result)) {
    PyObject *utf8 = PyUnicode_AsUTF8String(result);
    Py_DECREF(result);
    if (utf8) {
      Py_ssize_t len = PyBytes_GET_SIZE(utf8);
      memcpy(buffer, PyBytes_AS_STRING(utf8), MIN(len, bufsize));
      Py_DECREF(utf8);
      *ierr = 0;
      return (size_t)len;
    } else {
      *ierr = 1;
      return 0;
    }
  } else {
    Py_DECREF(result);
    PyErr_SetString(PyExc_TypeError, "Python file-like object read method "
                    "should return a 'str' object");
    *ierr = 1;
    return 0;
  }
}

static void pyfile_strreadlinefunc(void *data, GString *str, int *ierr)
{
  struct FileMethods *methods = (struct FileMethods *)data;
  static char fmt[] = "()";

  PyObject *result = PyObject_CallFunction(methods->readline_method, fmt);
  if (!result) {
    *ierr = 1;
  } else if (PyUnicode_Check(result)) {
#if PY_VERSION_HEX >= 0x03030000
    char *utf8 = PyUnicode_AsUTF8(result);
    if (utf8) {
      g_string_assign(str, utf8);
#else
    PyObject *utf8 = PyUnicode_AsUTF8String(result);
    if (utf8) {
      g_string_assign(str, PyBytes_AS_STRING(utf8));
      Py_DECREF(utf8);
#endif
      *ierr = 0;
    } else {
      *ierr = 1;
    }
    Py_DECREF(result);
  } else {
    Py_DECREF(result);
    PyErr_SetString(PyExc_TypeError, "Python file-like object readline method "
                    "should return a 'str' object");
    *ierr = 1;
  }
}

static void pyfile_strwritefunc(void *data, const char *buffer, size_t bufsize,
                                int *ierr)
{
  /* If an exception already occurred, exit immediately */
  if (PyErr_Occurred()) {
    *ierr = 1;
  } else {
    struct FileMethods *methods = (struct FileMethods *)data;
    static char fmt[] = "(s#)";

    PyObject *result = PyObject_CallFunction(methods->write_method, fmt,
                                             buffer,
#if PY_VERSION_HEX >= 0x02050000 && defined(PY_SSIZE_T_CLEAN)
                                             (Py_ssize_t)bufsize);
#else
                                             (int)bufsize);
#endif
    if (!result) {
      *ierr = 1;
    } else {
      *ierr = 0;
      Py_DECREF(result);
    }
  }
}
#endif

static void pyfile_freefunc(void *data)
{
  struct FileMethods *methods = (struct FileMethods *)data;
  Py_DECREF(methods->read_method);
  Py_DECREF(methods->readline_method);
  Py_DECREF(methods->write_method);
  g_free(methods);
}
%}

%inline %{
/* Open a mod_file object that accesses a Python file-like object */
static struct mod_file *mod_file_open_python(PyObject *obj, int *ierr)
{
  PyObject *read_method, *readline_method, *write_method;
  struct FileMethods *methods;

#if PY_VERSION_HEX < 0x03000000
  /* On Python 2 systems, try to use the faster cStringIO C API if available */
  if (PycStringIO && (PycStringIO_OutputCheck(obj)
                      || PycStringIO_InputCheck(obj))) {
    *ierr = 0;
    Py_INCREF(obj);
    return mod_file_open_stream(cstr_readfunc, cstr_readlinefunc,
                                cstr_writefunc, cstr_freefunc, obj);
  } else {
#endif
  read_method = PyObject_GetAttrString(obj, "read");
  if (!read_method) {
    PyErr_SetString(PyExc_TypeError,
                    "Expecting an object with a 'read' method");
    *ierr = 1;
    return NULL;
  }

  readline_method = PyObject_GetAttrString(obj, "readline");
  if (!readline_method) {
    Py_DECREF(read_method);
    PyErr_SetString(PyExc_TypeError,
                    "Expecting an object with a 'readline' method");
    *ierr = 1;
    return NULL;
  }

  write_method = PyObject_GetAttrString(obj, "write");
  if (!write_method) {
    Py_DECREF(read_method);
    Py_DECREF(readline_method);
    PyErr_SetString(PyExc_TypeError,
                    "Expecting an object with a 'write' method");
    *ierr = 1;
    return NULL;
  }

  methods = g_malloc(sizeof(struct FileMethods));
  methods->read_method = read_method;
  methods->readline_method = readline_method;
  methods->write_method = write_method;
  *ierr = 0;
#if PY_VERSION_HEX >= 0x03000000
  if (text_io_base && !PyObject_IsInstance(obj, text_io_base)) {
#endif
    return mod_file_open_stream(pyfile_readfunc, pyfile_readlinefunc,
                                pyfile_writefunc, pyfile_freefunc, methods);
#if PY_VERSION_HEX >= 0x03000000
  } else {
    return mod_file_open_stream(pyfile_strreadfunc, pyfile_strreadlinefunc,
                                pyfile_strwritefunc, pyfile_freefunc, methods);
  }
#endif
#if PY_VERSION_HEX < 0x03000000
  }
#endif
}
%} /* inline */

%include "../include/mod_file.h"
