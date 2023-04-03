import xdrlib
import os
import sys

try:
    import cPickle as pickle
except ImportError:
    import pickle


class TransferFile(object):
    """Object used to transfer a file between Communicators"""
    def __init__(self, filename):
        self._filename = filename


class HeartBeat:
    pass


class DataType:
    desc = "generic data"

    def recv(self, buffer):
        pass

    def send(self, obj):
        pass

    def __str__(self):
        return self.desc


if sys.version_info[0] >= 3:
    def _raw_byte_code(s):
        # Cannot simply use b'foo' since that is invalid syntax in Python 2.3
        return bytes(s, encoding='ascii')
else:
    def _raw_byte_code(s):
        return s

if sys.version_info[0] >= 3:
    # In Python 3, strings are Unicode, so they must be encoded as bytes
    # to travel across the network
    class NetString(DataType):
        code = _raw_byte_code('S')
        obj = str
        desc = "string"

        def recv(self, buffer):
            p = xdrlib.Unpacker(buffer[1:])
            obj = p.unpack_string().decode('UTF8')
            return (obj, buffer[1+p.get_position():])

        def send(self, obj):
            p = xdrlib.Packer()
            p.pack_string(obj.encode('UTF8'))
            return self.code + p.get_buffer()
else:
    class NetString(DataType):
        code = _raw_byte_code('S')
        obj = str
        desc = "string"

        def recv(self, buffer):
            p = xdrlib.Unpacker(buffer[1:])
            obj = p.unpack_string()
            return (obj, buffer[1+p.get_position():])

        def send(self, obj):
            p = xdrlib.Packer()
            p.pack_string(obj)
            return self.code + p.get_buffer()


class NetInteger(DataType):
    code = _raw_byte_code('I')
    obj = int
    desc = "integer"

    def recv(self, buffer):
        p = xdrlib.Unpacker(buffer[1:])
        obj = p.unpack_int()
        return (obj, buffer[1+p.get_position():])

    def send(self, obj):
        p = xdrlib.Packer()
        p.pack_int(obj)
        return self.code + p.get_buffer()


class NetFloat(DataType):
    code = _raw_byte_code('F')
    obj = float
    desc = "floating point"

    def recv(self, buffer):
        p = xdrlib.Unpacker(buffer[1:])
        obj = p.unpack_float()
        return (obj, buffer[1+p.get_position():])

    def send(self, obj):
        try:
            p = xdrlib.Packer()
            p.pack_float(obj)
            return self.code + p.get_buffer()
        # Python < 2.5 can fail trying to send Inf or NaN floats, so fall back
        # to pickling in this case:
        except SystemError:
            return NetPickle.send(obj)


class NetPickle(DataType):
    code = _raw_byte_code('P')
    obj = None

    desc = "Python pickled object"

    def recv(self, buffer):
        p = xdrlib.Unpacker(buffer[1:])
        obj = p.unpack_string()
        return (pickle.loads(obj), buffer[1+p.get_position():])

    def send(self, obj):
        p = xdrlib.Packer()
        try:
            p.pack_string(pickle.dumps(obj, -1))
        # Python < 2.5 can fail trying to send Inf or NaN floats in binary
        # mode, so fall back to the old protocol in this case:
        except SystemError:
            p.pack_string(pickle.dumps(obj, 0))
        return self.code + p.get_buffer()


class NetFile(DataType):
    code = _raw_byte_code('f')
    obj = TransferFile
    desc = "Transferred file"

    def recv(self, buffer):
        filename, buffer = NetString.recv(buffer[1:])
        filelen, buffer = NetInteger.recv(buffer)
        if len(buffer) < filelen:
            raise IndexError("File not completely read")
        filename = os.path.basename(filename)
        f = open(filename, 'wb')
        f.write(buffer[:filelen])
        return TransferFile(filename), buffer[filelen:]

    def send(self, obj):
        filename = obj._filename
        f = open(filename, 'rb')
        buf = f.read()
        return (self.code + NetString.send(filename)
                + NetInteger.send(len(buf)) + buf)


class NetCommand(NetString):
    code = _raw_byte_code('C')
    desc = "command"


typemap = {}
cmdmap = {}
for name in dir():
    var = eval(name)
    try:
        if issubclass(var, DataType) and var is not DataType:
            exec("%s = %s()" % (name, name))
            var = eval(name)
            if sys.version_info[0] >= 3:
                # In Python 3, b'foo'[0] returns an int, not b'f'
                cmdmap[ord(var.code)] = var
            else:
                cmdmap[var.code] = var
            typemap[var.obj] = var
    except TypeError:
        pass
