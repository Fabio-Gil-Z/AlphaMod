from modeller.parallel.worker import Worker
from modeller.parallel.myspawn import myspawn
from modeller.util.deprecation import _deprecation_handler


class SSHWorker(Worker):
    """A worker on a remote host accessed via ssh or rsh"""

    def __init__(self, nodename, ssh_command='ssh'):
        self._nodename = nodename
        self._ssh = ssh_command
        Worker.__init__(self)

    def start(self, path, id, output):
        Worker.start(self, path, id, output)
        myspawn("%s %s %s -worker %s" % (self._ssh, self._nodename, path, id),
                output)

    def __repr__(self):
        return "<ssh worker on %s>" % self._nodename


# Compatibility with Modeller 9
class ssh_slave(SSHWorker):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(ssh_slave)
        SSHWorker.__init__(self, *args, **keys)
