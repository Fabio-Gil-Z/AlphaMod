from modeller.parallel.worker import Worker
from modeller.parallel.myspawn import myspawn
from modeller.util.deprecation import _deprecation_handler


class LocalWorker(Worker):
    """A worker running on the same machine"""

    def start(self, path, id, output):
        Worker.start(self, path, id, output)
        myspawn("%s -worker %s" % (path, id), output)

    def __repr__(self):
        return "<Worker on localhost>"


# Compatibility with Modeller 9
class local_slave(LocalWorker):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(local_slave)
        LocalWorker.__init__(self, *args, **keys)
