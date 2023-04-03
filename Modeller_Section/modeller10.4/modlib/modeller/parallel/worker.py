import time
from modeller.parallel import communicator, data_types, workerstate
from modeller.util.deprecation import _deprecation_handler


class Worker(communicator.Communicator):

    def __init__(self):
        communicator.Communicator.__init__(self)
        self.task = None
        self._state = workerstate.init
        self.update_contact_time()

    def start(self, path, id, output):
        """Start the worker running on the remote host"""
        self._state = workerstate.pending

    def accept_connection(self, socket):
        communicator.Communicator.accept_connection(self, socket)
        self._state = workerstate.connected
        self.update_contact_time()

    def run_cmd(self, cmd):
        """Run a Python command on the remote host"""
        self._send(data_types.NetCommand.send(cmd))

    def set_directory(self, dir):
        self.run_cmd('import os\nos.chdir(manager.get_data())')
        self.send_data(dir)

    def set_python_search_path(self, path):
        self.run_cmd('import sys\nsys.path.insert(0, manager.get_data())')
        self.send_data(path)

    def set_log_level(self, log):
        for typ in ('output', 'notes', 'warnings', 'errors', 'memory'):
            self.run_cmd('log.%s = manager.get_data()' % typ)
            self.send_data(eval("log.%s" % typ))

    def update_contact_time(self):
        self.last_contact_time = time.time()

    def contact_timeout(self, timeout):
        return (time.time() - self.last_contact_time) > timeout

    def run_task(self, task):
        if not self.ready_for_task():
            raise TypeError("%s not ready for task" % str(self))
        self._state = workerstate.running_task
        self.task = task
        self.run_cmd('task = manager.get_data()')
        self.run_cmd('task._setup()')
        self.send_data(task)
        for transfer in task.input_files:
            self.send_data(data_types.TransferFile(transfer))
        self.run_cmd('task._results = task._do_run(manager)')
        self.run_cmd('manager.send_data(task._results)')
        self.run_cmd('del task')

    def task_results(self):
        """Read results from worker. If these terminate the task, return it."""
        while True:
            r = self.get_data(allow_heartbeat=True)
            self.update_contact_time()
            # Allow for the case where both a heartbeat and task results are in
            # the receive buffer:
            if isinstance(r, data_types.HeartBeat):
                if not self.data_pending():
                    return None
            else:
                break
        task = self.task
        task._results = r
        self.task = None
        self._state = workerstate.connected
        return task

    def kill(self):
        self.disconnect()
        task = self.task
        self.task = None
        self._state = workerstate.dead
        return task

    def get_state(self):
        return self._state

    def ready_for_task(self):
        return self.get_state() == workerstate.connected

    def running_task(self):
        return self.get_state() == workerstate.running_task


# Compatibility with Modeller 9
class slave(Worker):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(slave)
        Worker.__init__(self, *args, **keys)
