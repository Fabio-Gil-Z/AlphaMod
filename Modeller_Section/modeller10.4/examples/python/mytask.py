from modeller import *
from modeller.parallel import Task

class MyTask(Task):
    """A task to read in a PDB file on the worker, and return the resolution"""
    def run(self, code):
        env = Environ()
        env.io.atom_files_directory = ["../atom_files"]
        mdl = Model(env, file=code)
        return mdl.resolution
