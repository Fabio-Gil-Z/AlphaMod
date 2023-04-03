class WorkerState:
    def __init__(self, desc):
        self.desc = desc


init = WorkerState("Initial state")
pending = WorkerState("Pending startup")
connected = WorkerState("Connected, with no active task")
running_task = WorkerState("Connected, and running a task")
dead = WorkerState("Disconnected due to failure")
