from modeller.parallel.job import Job
from modeller.parallel.sge_qsub_array import SGEQsubArray
from modeller.util.deprecation import _deprecation_handler


class SGEQsubJob(Job):
    """A parallel job which automatically starts Sun Grid Engine processes
       using 'qsub'."""

    def __init__(self, options, maxworker, seq=(), modeller_path=None,
                 host=None):
        Job.__init__(self, seq, modeller_path, host)
        self.options = options
        self.maxworker = maxworker
        self.pending_arrays = []

    def expand_for_tasks(self):
        numworkers = len(self)
        numdesired = len(self.tasks)
        if self.maxworker is not None and self.maxworker < numdesired:
            numdesired = self.maxworker
        if numworkers < numdesired:
            array = SGEQsubArray(self.options, numdesired - numworkers)
            self.extend(array)
            self.pending_arrays.append(array)

    def start_processes(self, port):
        """Start all new workers"""
        Job.start_processes(self, port)
        jobname = self.get_name()
        for array in self.pending_arrays:
            array.start(jobname)
        self.pending_arrays = []


# Compatibility with Modeller 9
class sge_qsub_job(SGEQsubJob):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(sge_qsub_job)
        SGEQsubJob.__init__(self, *args, **keys)
