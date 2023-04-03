import socket
import sys
from modeller.parallel import communicator, data_types
from modeller.parallel.worker_communicator import WorkerCommunicator
import threading


class HeartBeat(threading.Thread):
    """Periodically send a 'heartbeat' back to the manager, so that it can
       distinguish between failed nodes and long calculations"""
    timeout = 300

    def __init__(self, manager):
        threading.Thread.__init__(self)
        self.manager = manager
        self.event = threading.Event()

    def cancel(self):
        """Stop the heartbeat"""
        self.event.set()

    def run(self):
        while True:
            self.event.wait(self.timeout)
            if self.event.isSet():
                break
            else:
                self.manager.send_data(data_types.HeartBeat())


class WorkerLoop(object):
    def __init__(self, addr):
        self.addr = addr

    def _handle_worker_io(self, manager, workerdict):
        """Handle all messages from manager"""
        while True:
            try:
                cmdstr = manager.get_command()
            except communicator.NetworkError:
                # Connection broken - shutdown worker
                break
            try:
                exec(cmdstr, workerdict)
            except Exception:
                detail = sys.exc_info()[1]
                # Propagate errors to manager, and reraise them here (but don't
                # send back erorrs we got from the manager!)
                if not isinstance(detail, communicator.RemoteError):
                    try:
                        manager.send_data(communicator.ErrorWrapper(detail))
                    except socket.error:
                        detail2 = sys.exc_info()[1]
                        print("Warning: ignored exception " + str(detail2)
                              + " when trying to send error state "
                              + str(detail) + " back to manager")
                        raise detail
                raise

    def run(self):
        print("Worker startup: connect to manager at %s" % self.addr)
        (host, port, identifier) = self.addr.split(":", 2)
        port = int(port)
        lock = threading.Lock()
        manager = WorkerCommunicator(lock, reconnect=(host, port, identifier))
        manager.connect_back(host, port, identifier)
        th = HeartBeat(manager)
        th.start()
        # 'master' is for Modeller 9 compatibility
        workerdict = {'manager': manager, 'master': manager}
        exec('from modeller import *', workerdict)
        try:
            self._handle_worker_io(manager, workerdict)
        finally:
            th.cancel()
        print("Worker shutdown")
