"""Classes used by AutoModel and LoopModel to parallelize model building"""

from modeller.parallel.task import Task


class ModelTask(Task):
    """Build a single model in parallel"""
    def __init__(self, automdl, num, atmsel):
        Task.__init__(self, automdl, num, atmsel)
        self.num = num

    def __str__(self):
        return "<Model building task #%d>" % self.num

    def run(self, automdl, num, atmsel):
        # Make sure we get a different initial structure for each model
        automdl.env.libs.random_perturb(num)
        return automdl.single_model(atmsel, num, True)


class LoopTask(Task):
    """Build a single loop model in parallel"""
    def __init__(self, loopmdl, atmsel, ini_model, num, id1, sched):
        Task.__init__(self, loopmdl, atmsel, ini_model, num, id1, sched)
        self.num = num
        self.id1 = id1

    def __str__(self):
        return "<Loop model building task #%d.%d>" % (self.num, self.id1)

    def run(self, loopmdl, atmsel, ini_model, num, id1, sched):
        # Make sure we get a different initial structure for each model
        loopmdl.env.libs.random_perturb(id1)
        return loopmdl.single_loop_model(atmsel, ini_model, num, id1, sched,
                                         True)
