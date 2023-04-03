"""A simple script to start up a Modeller worker process. This is useful if,
   for example, you want to use your system's native Python version rather
   than the Python 2.3 interpreter compiled into Modeller. In most cases this
   is handled automatically.

   Use by setting the 'modeller_path' variable when creating a 'job' object,
   to something like
   "/usr/bin/python <modeller_dir>/modlib/modeller/parallel/modworker.py".
   (Note that if your Modeller Python modules are not
   in the default search path, you may have to write a wrapper script to set
   PYTHONPATH and/or LD_LIBRARY_PATH or similar so that the modules can be
   located.)
"""

import sys


def usage():
    print("Usage: %s -worker managerspec\n" % sys.argv[0])
    sys.exit(1)


def main():
    if len(sys.argv) != 3 or sys.argv[1] != '-worker':
        usage()
    else:
        from modeller.parallel.workerloop import WorkerLoop
        # Ensure that the current directory is in the search path:
        sys.path.insert(0, '')
        s = WorkerLoop(sys.argv[2])
        s.run()


if __name__ == '__main__':
    main()
