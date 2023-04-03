import sys
try:
    import subprocess
except ImportError:
    from modeller.python_library import subprocess


def myspawn(cmd, output):
    """Run ``cmd`` in the background, and direct stdout and stderr to
       ``output``."""

    fp = open(output, "w")
    print("%s >& %s" % (cmd, output))
    if sys.platform == 'win32':
        _myspawn_win32(cmd, fp)
    else:
        _myspawn_unix(cmd, fp)


def _myspawn_unix(cmd, fp):
    p = subprocess.Popen(cmd, shell=True, stdout=fp, stderr=subprocess.STDOUT)


def _myspawn_win32(cmd, fp):
    try:
        # shell isn't needed on Win32, and may not be found under wine anyway
        p = subprocess.Popen(cmd, shell=False, stdout=fp,
                             stderr=subprocess.STDOUT)
    # Ignore Windows "file not found" errors, so that behavior is consistent
    # between Unix and Windows
    except WindowsError:
        print("WindowsError: %s (ignored)" % sys.exc_info()[1])
