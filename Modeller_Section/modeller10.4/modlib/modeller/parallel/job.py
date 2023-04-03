import socket
import random
import os
import select
import errno
import sys
import modeller
import _modeller
from modeller.parallel.communicator import NetworkError, TaskSetupError
from modeller.parallel import workerstate
from modeller.util.deprecation import _deprecation_handler


def _ignore_eintr(func, *args, **keys):
    "Call the given function, and ignore any 'interrupted system call' error"
    while True:
        try:
            return func(*args, **keys)
        except socket.timeout:
            raise  # pass timeout exceptions unchanged
        except (select.error, socket.error):
            exc = sys.exc_info()[1]
            if exc[0] != errno.EINTR:
                raise


class Job(list):
    connect_timeout = 7200
    heartbeat_timeout = 7200

    def __init__(self, seq=(), modeller_path=None, host=None):
        list.__init__(self, seq)
        self.worker_startup_commands = []
        self.tasks = []
        if modeller_path is not None:
            self.modeller_path = modeller_path
        else:
            self.modeller_path = self.get_default_modeller_path()
        if host:
            self.host = host
        else:
            self.host = self.__get_default_hostname()
        self.listensock = self.__listen_on_random_port(self.host,
                                                       self.connect_timeout)
        self.pending_workers = {}
        self.connected_workers = {}
        self.cwd = os.getcwd()

    # Modeller 9 compatibility
    slave_startup_commands = property(
        lambda self: self.worker_startup_commands)

    def __get_default_hostname(self):
        try:
            # Get primary IP address of this machine
            return socket.gethostbyname_ex(socket.gethostname())[-1][0]
        except socket.gaierror:
            print("""
Could not determine hostname (this usually indicates that networking is not
set up correctly). Defaulting to 'localhost'. Note that this will prevent any
worker except for LocalWorker from working correctly. If you need to use such
workers, set the 'host' parameter manually when creating the Job object.
""")
        return 'localhost'

    def get_default_modeller_path(self):
        def quote(pth):
            # quote the path if it contains spaces so as to not confuse
            # the shell
            if ' ' in pth:
                return '"%s"' % pth
            else:
                return pth
        bindir = modeller.info.bindir
        modpy = os.path.join(bindir, 'modpy.sh')
        try:
            from modeller.parallel import modworker
        except ImportError:
            modworker = None
        if modworker and not _modeller.mod_embedded_get():
            path = quote(sys.executable) + " " + quote(modworker.__file__)
            if os.path.exists(modpy):
                path = modpy + " " + path
            return path
        else:
            return "mod" + _modeller.mod_short_version_get()

    def get_name(self):
        job = _modeller.mod_jobname_get()
        if job == '(stdin)':
            job = 'stdout'
        else:
            job = os.path.basename(job)
        return job

    def start_processes(self, port):
        job = self.get_name()
        for (num, worker) in enumerate(self):
            if worker.get_state() == workerstate.init:
                id = self.__get_id(num)
                self.pending_workers[id] = worker
                addr = "%s:%d" % (self.host, port)
                output = job + ".worker%d" % num
                worker.start(self.modeller_path, "%s:%s" % (addr, id), output)

    def expand_for_tasks(self):
        pass

    def accept_worker(self, sock, pending_workers, connected_workers):
        # Make sure the new socket is blocking (on some platforms this socket
        # inherits non-blocking status from the listening socket)
        sock.setblocking(True)
        id = sock.recv(1024)
        if sys.version_info[0] >= 3:
            id = id.decode('ascii')
        if id and id in pending_workers:
            worker = pending_workers.pop(id)
            connected_workers[id] = worker
            print("Identified worker %s " % str(worker))
            worker.accept_connection(sock)
            worker.set_directory(self.cwd)
            if sys.path[0] != '':
                worker.set_python_search_path(sys.path[0])
            for cmd in self.worker_startup_commands:
                worker.run_cmd(cmd)
            worker.set_log_level(modeller.log)
            return worker
        elif id and id in connected_workers:
            worker = connected_workers[id]
            print("Reconnect from worker %s " % str(worker))
            worker.accept_connection(sock)
        else:
            print("Ignoring request from unknown worker")

    def start(self):
        """Start all non-running workers"""
        (s, port) = self.listensock
        self.start_processes(port)
        while len(self.pending_workers) > 0:
            (conn, addr) = _ignore_eintr(s.accept)
            self.accept_worker(conn, self.pending_workers,
                               self.connected_workers)
        print("All workers connected OK")

    def __listen_on_random_port(self, host, timeout):
        """Open a listening socket on a random high-numbered port"""
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        tries = 0
        while True:
            port = random.randint(10000, 60000)
            try:
                s.bind((host, port))
            # gaierror is a subclass of error, so catch it separately
            except socket.gaierror:
                raise
            except socket.error:
                tries += 1
                if tries > 10:
                    raise
            else:
                break
        s.listen(15)
        s.settimeout(timeout)
        return (s, port)

    def queue_task(self, taskobj):
        self.tasks.append(taskobj)

    def run_all_tasks(self):
        """Run all tasks and return all the results, in the same order that they
           were submitted, when all the jobs have completed."""
        tasks = self.tasks[:]
        print("Running %d tasks on %d workers" % (len(tasks), len(self)))
        self.push_tasks_to_workers()
        while True:
            try:
                for task in self._finish_all_tasks():
                    pass
            except IndexError:
                break
        if len(self.tasks) > 0:
            raise ValueError("Ran out of workers to run tasks")
        return [task._results for task in tasks]

    def yield_tasks_unordered(self):
        """Run all tasks and return their results (as a generator), in
           whatever order they complete."""
        print("Running %d tasks on %d workers" % (len(self.tasks), len(self)))
        self.push_tasks_to_workers()
        while True:
            try:
                for task in self._finish_all_tasks():
                    yield task._results
            except IndexError:
                break
        if len(self.tasks) > 0:
            raise ValueError("Ran out of workers to run tasks")

    def _finish_all_tasks(self):
        """Waits for tasks to finish from any running worker; if any pending
           workers try to connect, start them up. Return each finished task,
           as a generator."""
        (s, port) = self.listensock
        while True:
            events = self.get_next_events()
            if events is None:
                self.kill_all_running_workers()
            else:
                for obj in events:
                    task = self._process_event(obj, s)
                    if task:
                        yield task

    def _process_event(self, obj, listensock):
        """Handle a single event returned from a worker"""
        if obj == listensock:
            # New worker just connected to the listening socket
            (conn, addr) = _ignore_eintr(listensock.accept)
            worker = self.accept_worker(conn, self.pending_workers,
                                        self.connected_workers)
            if worker and len(self.tasks) > 0:
                worker.run_task(self.tasks.pop(0))
        elif obj.running_task():
            # A worker returned data
            try:
                task = obj.task_results()
                if task:
                    # The worker completed its task
                    print("%s on %s completed" % (str(task), str(obj)))
                    if len(self.tasks) > 0:
                        obj.run_task(self.tasks.pop(0))
                    return task
                else:
                    # The worker sent back a heartbeat; check for any
                    # dead workers
                    self.kill_timed_out_workers()
            except (NetworkError, TaskSetupError):
                self.kill_workers((obj,), sys.exc_info()[1])
        else:
            print("Warning: worker %s reports data, but is not running a task"
                  % str(obj))

    def kill_workers(self, workers, err=""):
        if err != "":
            err = "(%s) " % err
        for s in workers:
            print("%s failed %s- removing from %s" % (s, err, self))
            task = s.kill()
            if task:
                self.tasks.append(task)
        self.push_tasks_to_workers()

    def kill_all_running_workers(self):
        running = [a for a in self if a.running_task()]
        self.kill_workers(running)
        raise NetworkError("Did not hear from any running worker in %d seconds"
                           % self.heartbeat_timeout)

    def kill_timed_out_workers(self):
        timedout = [a for a in self if a.running_task() and
                    a.contact_timeout(self.heartbeat_timeout)]
        if len(timedout) > 0:
            print("Did not hear from workers %s in %d seconds" %
                  (str(timedout), self.heartbeat_timeout))
            self.kill_workers(timedout)

    def push_tasks_to_workers(self):
        (s, port) = self.listensock
        self.start_processes(port)
        for worker in [a for a in self if a.ready_for_task()]:
            try:
                t = self.tasks.pop(0)
            except IndexError:
                break
            try:
                worker.run_task(t)
            # If a network error occurred, kill the worker and requeue the task
            except socket.error:
                print("worker %s failed on run task with %s; "
                      "removing from job" % (worker, sys.exc_info()[1]))
                worker.kill()
                self.tasks.insert(0, t)
        self.expand_for_tasks()
        self.start_processes(port)

    def get_next_events(self):
        (s, port) = self.listensock
        running = [a for a in self if a.running_task()]
        workermap = {}
        fileno = s.fileno()
        workermap[fileno] = s
        if len(running) == 0:
            if len(self.tasks) == 0:
                raise IndexError("No more tasks")
            if len(self.pending_workers) == 0:
                raise IndexError("No more workers")
        # poll() for new events, or fall back to select() on platforms
        # which don't have poll().
        try:
            poll = select.poll()
        except AttributeError:
            poll = None
        if poll:
            return self.__next_events_poll(poll, running, workermap, fileno)
        else:
            return self.__next_events_select(running, workermap, fileno)

    def __next_events_poll(self, poll, running, workermap, fileno):
        poll.register(fileno, select.POLLIN)
        for worker in running:
            fileno = worker.socket.fileno()
            workermap[fileno] = worker
            poll.register(fileno, select.POLLIN)
        ready = poll.poll(self.heartbeat_timeout * 1000)
        if len(ready) == 0:
            return None
        else:
            return [workermap[fd[0]] for fd in ready]

    def __next_events_select(self, running, workermap, fileno):
        waitin = [fileno]
        for worker in running:
            fileno = worker.socket.fileno()
            workermap[fileno] = worker
            waitin.append(fileno)
        (ready, rout, rerr) = _ignore_eintr(select.select, waitin, [], [],
                                            self.heartbeat_timeout)
        if len(ready) == 0:
            return None
        else:
            return [workermap[fd] for fd in ready]

    def __get_id(self, num):
        """Return a random identifier, used to make sure the right workers
           connect back to us."""
        id = "%d:" % num
        for i in range(0, 8):
            id += chr(random.randint(0, 25) + ord('A'))
        return id

    def __repr__(self):
        return "<Parallel job [" + \
               ", ".join([str(node) for node in self]) + "]>"


# Compatibility with Modeller 9
class job(Job):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(job)
        Job.__init__(self, *args, **keys)
