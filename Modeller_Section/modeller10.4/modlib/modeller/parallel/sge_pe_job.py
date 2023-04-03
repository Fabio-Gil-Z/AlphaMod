import os
from modeller.parallel.job import Job
from modeller.parallel.sge_pe_worker import SGEPEWorker
from modeller.parallel.local_worker import LocalWorker
from modeller.util.deprecation import _deprecation_handler


class SGEPEJob(Job):
    """A parallel job containing processes on all Sun Grid Engine
       worker nodes"""

    def __init__(self, seq=(), modeller_path=None, host=None):
        Job.__init__(self, seq, modeller_path, host)
        pe = os.environ['PE_HOSTFILE']
        fh = open(pe, "r")
        while True:
            line = fh.readline()
            if line == '':
                break
            (node, num, queue) = line.split(None, 2)
            for i in range(int(num)):
                self.append(SGEPEWorker(node))
        # Replace first worker with a LocalWorker, as this is ourself, and SGE
        # won't let us start this process with qrsh (as we are already
        # occupying the slot)
        if len(self) > 0:
            self[0] = LocalWorker()


# Compatibility with Modeller 9
class sge_pe_job(SGEPEJob):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(sge_pe_job)
        SGEPEJob.__init__(self, *args, **keys)
