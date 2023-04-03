from modeller.parallel import data_types
import os
import tempfile
import shutil
import sys
from modeller.parallel.communicator import TaskSetupError
from modeller.util.deprecation import _deprecation_handler


class Task(object):
    _results = None
    run_in_tempdir = False

    def __init__(self, *args, **vars):
        self.input_files = []
        self.output_files = []
        self._args = args
        self._vars = vars

    def __str__(self):
        return "<Task>"

    def _setup(self):
        """Do pre-run setup, e.g. making a temporary run directory"""
        try:
            if self.run_in_tempdir:
                self._cwd = os.getcwd()
                self._tmpdir = tempfile.mkdtemp()
                os.chdir(self._tmpdir)
        except Exception:
            # Wrap errors that occur in the setup phase
            raise TaskSetupError(sys.exc_info()[1])

    def _do_run(self, manager):
        """Actually run the task"""
        try:
            ret = self.run(*self._args, **self._vars)
            for transfer in self.output_files:
                manager.send_data(data_types.TransferFile(transfer))
        finally:
            if self.run_in_tempdir:
                os.chdir(self._cwd)
                shutil.rmtree(self._tmpdir, ignore_errors=True)
        return ret

    def run(self):
        raise NotImplementedError


# Compatibility with Modeller 9
class task(Task):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(task)
        Task.__init__(self, *args, **keys)
