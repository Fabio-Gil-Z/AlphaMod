from modeller.parallel.communicator import Communicator
import socket
import sys


class WorkerCommunicator(Communicator):
    connect_timeout = 600

    def __init__(self, lock=None, reconnect=None):
        Communicator.__init__(self, lock)
        self.reconnect = reconnect

    def connect_back(self, host, port, identifier):
        """Establish a TCP connection with the manager"""
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(self.connect_timeout)
        s.connect((host, port))
        if sys.version_info[0] >= 3:
            # Strings are Unicode in Python 3; must encode for transport
            s.sendall(identifier.encode('ascii'))
        else:
            s.sendall(identifier)
        s.settimeout(None)
        self.accept_connection(s)

    def _send(self, data):
        """Try to reopen the connection to the manager if we got a broken
           pipe"""
        try:
            self.socket.sendall(data)
        except socket.error:
            if self.reconnect:
                self.connect_back(*self.reconnect)
                self.socket.sendall(data)
            else:
                raise
