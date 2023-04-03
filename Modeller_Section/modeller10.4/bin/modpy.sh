#!/bin/sh

TOPDIR="/home/pc7/bin/modeller10.4"
EXETYPE=x86_64-intel8

LLP=${TOPDIR}/lib/${EXETYPE}
if test -z "${LD_LIBRARY_PATH}"; then
  LD_LIBRARY_PATH=${LLP}
else
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${LLP}
fi
if test -z "${DYLD_LIBRARY_PATH}"; then
  DYLD_LIBRARY_PATH=${LLP}
else
  DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${LLP}
fi
if test -z "${LIBPATH}"; then
  LIBPATH=${LLP}
else
  LIBPATH=${LIBPATH}:${LLP}
fi

if test "x${EXETYPE}" = "xi386-w32" -o "x${EXETYPE}" = "xx86_64-w64"; then
  PP=${TOPDIR}/modlib
else
  PP=${TOPDIR}/lib/${EXETYPE}:${TOPDIR}/modlib
fi
if test -z "${PYTHONPATH}"; then
  PYTHONPATH=${PP}
else
  ORIGPYPATH="${PYTHONPATH}"
  if test "x${EXETYPE}" = "xi386-w32" -o "x${EXETYPE}" = "xx86_64-w64"; then
    PYTHONPATH="${PYTHONPATH};${PP}"
  else
    PYTHONPATH=${PYTHONPATH}:${PP}
  fi
fi
export LD_LIBRARY_PATH DYLD_LIBRARY_PATH PYTHONPATH LIBPATH ORIGPYPATH
exec "$@"
