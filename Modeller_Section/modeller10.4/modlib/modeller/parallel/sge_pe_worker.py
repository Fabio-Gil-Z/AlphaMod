from modeller.parallel.worker import Worker
from modeller.parallel.myspawn import myspawn
from modeller.util.deprecation import _deprecation_handler


class SGEPEWorker(Worker):
    """A worker within a Sun Grid Engine parallel environment"""

    def __init__(self, nodename):
        self._nodename = nodename
        Worker.__init__(self)

    def start(self, path, id, output):
        Worker.start(self, path, id, output)
        myspawn("qrsh -inherit -V %s %s -worker %s"
                % (self._nodename, path, id), output)

    def __repr__(self):
        return "<SGE PE worker on %s>" % self._nodename


# Compatibility with Modeller 9
class sge_pe_slave(SGEPEWorker):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(sge_pe_slave)
        SGEPEWorker.__init__(self, *args, **keys)
