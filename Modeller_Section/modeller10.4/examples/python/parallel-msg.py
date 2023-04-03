from modeller import *
from modeller.parallel import Job, LocalWorker

# Create an empty parallel job, and then add a single worker process running
# on the local machine
j = Job()
j.append(LocalWorker())

# Start all worker processes (note: this will only work if 'modxxx' - where
# xxx is the Modeller version - is in the PATH; if not, use modeller_path
# to specify an alternate location)
j.start()

# Have each worker read in a PDB file (provided by us, the master) and
# return the PDB resolution back to us
for worker in j:
    worker.run_cmd('''
env = Environ()
env.io.atom_files_directory = ["../atom_files"]
log.verbose()
code = master.get_data()
mdl = Model(env, file=code)
master.send_data(mdl.resolution)
''')
    worker.send_data('1fdn')
    data = worker.get_data()
    print("%s returned model resolution: %f" % (str(worker), data))
