from modeller.parallel import data_types
import socket
import sys


class RemoteError(Exception):
    def __init__(self, exc, worker):
        self.exc = exc
        self.worker = worker

    def __str__(self):
        errstr = str(self.exc.__class__).replace("exceptions.", "")
        return "%s: %s from %s" % (errstr, str(self.exc), str(self.worker))


class TaskSetupError(Exception):
    pass


class ErrorWrapper(object):
    def __init__(self, exc):
        self.exc = exc


class NetworkError(Exception):
    pass


class Communicator(object):
    def __init__(self, lock=None):
        # Buffer must be raw bytes, not a Unicode string
        if sys.version_info[0] >= 3:
            self.buffer = bytes()
        else:
            self.buffer = ''
        self.socket = None
        self.lock = lock
        self.__pending_commands = []

    def _send(self, data):
        self.socket.sendall(data)

    def accept_connection(self, socket):
        """Use the given network socket to talk to the remote host"""
        self.socket = socket

    def disconnect(self):
        """Shut down the connection to the remote host"""
        self.acquire_lock()
        try:
            if self.socket:
                self.socket.close()
            self.socket = None
        finally:
            self.release_lock()

    def acquire_lock(self):
        """Acquire the lock to allow only one thread to access the network
           at a time"""
        if self.lock:
            self.lock.acquire()

    def release_lock(self):
        """Release the lock to allow only one thread to access the network
           at a time"""
        if self.lock:
            self.lock.release()

    def send_data(self, data):
        """Send some data to the remote host"""
        self.acquire_lock()
        try:
            try:
                obj = data_types.typemap[data.__class__]
            except KeyError:
                obj = data_types.NetPickle
            self._send(obj.send(data))
        finally:
            self.release_lock()

    def data_pending(self):
        """Return True if data has been read from the network and is available
           for get_data"""
        return len(self.buffer) > 0

    def get_data(self, allow_heartbeat=False):
        """Get some data from the remote host. If allow_heartbeat is False,
           any heartbeat data is ignored, and the next data packet is returned
           instead."""
        self.acquire_lock()
        try:
            while True:
                (cmdtype, obj) = self._recv()
                if cmdtype == data_types.NetCommand:
                    self.__pending_commands.append(obj)
                elif (isinstance(obj, data_types.HeartBeat)
                      and not allow_heartbeat):
                    pass
                else:
                    break
            return obj
        finally:
            self.release_lock()

    def get_command(self):
        self.acquire_lock()
        try:
            return self.__get_command_locked()
        finally:
            self.release_lock()

    def __get_command_locked(self):
        pending = self.__pending_commands
        if len(pending) > 0:
            return pending.pop(0)
        else:
            (cmdtype, obj) = self._recv()
            if cmdtype == data_types.NetCommand:
                return obj
            else:
                raise TypeError("Was expecting a command; got " + str(obj) +
                                " data instead")

    def _recv(self):
        while True:
            try:
                cmd = self.buffer[0]
                cmdtype = data_types.cmdmap[cmd]
                (obj, self.buffer) = cmdtype.recv(self.buffer)
                if isinstance(obj, ErrorWrapper):
                    if isinstance(obj.exc, TaskSetupError):
                        raise obj.exc
                    else:
                        raise RemoteError(obj.exc, self)
                elif isinstance(obj, data_types.TransferFile):
                    del obj  # file is already dumped automatically to disk
                else:
                    return (cmdtype, obj)
            except (IndexError, EOFError):
                try:
                    data = self.socket.recv(1024)
                except socket.error:
                    raise NetworkError("Connection lost to worker %s: %s"
                                       % (str(self), str(sys.exc_info()[1])))
                if len(data) == 0:
                    raise NetworkError("Connection lost to worker %s"
                                       % str(self))
                self.buffer += data
