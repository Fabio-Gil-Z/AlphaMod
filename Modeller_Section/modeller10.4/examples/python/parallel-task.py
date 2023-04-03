from modeller import *
from modeller.parallel import Job, LocalWorker

# Load in my task from mytask.py (note: needs to be in a separate Python
# module like this, in order for Python's pickle module to work correctly)
from mytask import MyTask

log.minimal()
# Create an empty parallel job, and then add 2 worker processes running
# on the local machine
j = Job()
j.append(LocalWorker())
j.append(LocalWorker())

# Run 'mytask' tasks
j.queue_task(MyTask('1fdn'))
j.queue_task(MyTask('1b3q'))
j.queue_task(MyTask('1blu'))

results = j.run_all_tasks()

print("Got model resolution: " + str(results))
