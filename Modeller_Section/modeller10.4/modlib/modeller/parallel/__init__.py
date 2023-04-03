from modeller.parallel.job import Job
from modeller.parallel.sge_pe_job import SGEPEJob
from modeller.parallel.sge_qsub_job import SGEQsubJob
from modeller.parallel.worker import Worker
from modeller.parallel.sge_pe_worker import SGEPEWorker
from modeller.parallel.sge_qsub_worker import SGEQsubWorker
from modeller.parallel.ssh_worker import SSHWorker
from modeller.parallel.local_worker import LocalWorker
from modeller.parallel.communicator import NetworkError, RemoteError
from modeller.parallel.task import Task

# Compatibility with Modeller 9
from modeller.parallel.job import job
from modeller.parallel.sge_pe_job import sge_pe_job
from modeller.parallel.sge_qsub_job import sge_qsub_job
from modeller.parallel.worker import slave
from modeller.parallel.sge_pe_worker import sge_pe_slave
from modeller.parallel.sge_qsub_worker import sge_qsub_slave
from modeller.parallel.ssh_worker import ssh_slave
from modeller.parallel.local_worker import local_slave
from modeller.parallel.task import task
