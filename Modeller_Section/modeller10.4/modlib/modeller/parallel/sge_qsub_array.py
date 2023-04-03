import re
from modeller.parallel.sge_qsub_worker import SGEQsubWorker
from modeller.parallel.mypopen4 import MyPopen4


class SGEQsubArray(list):
    """An array of workers started with Sun Grid Engine 'qsub'"""

    def __init__(self, options, numworker, seq=()):
        list.__init__(self, seq)
        self._options = options
        self._pending_workers = []
        for i in range(numworker):
            self.append(SGEQsubWorker(options, array=self))

    def start_worker(self, path, id, output, standard_options):
        self._path = path
        self._standard_options = standard_options
        self._pending_workers.append((id, output))

    def start(self, jobname):
        qsub = "qsub -S /bin/sh -N '%s' -o sge-errors %s %s -t 1-%d" % \
               (jobname, self._options, self._standard_options,
                len(self._pending_workers))
        print(qsub)
        a = MyPopen4(qsub)
        (input, output) = (a.stdin, a.stdout)
        id = " ".join([repr(s[0]) for s in self._pending_workers])
        out = " ".join([repr(s[1]) for s in self._pending_workers])
        input.write("#!/bin/sh\n")
        input.write("id=( '' %s )\n" % id)
        input.write("out=( '' %s )\n" % out)
        input.write("myid=${id[$SGE_TASK_ID]}\n")
        input.write("myout=${out[$SGE_TASK_ID]}\n")
        input.write("%s -worker $myid > $myout 2>&1\n" % self._path)
        input.close()
        outlines = output.readlines()
        output.close()
        for line in outlines:
            print(line.rstrip('\r\n'))
        a.require_clean_exit()
        self._set_jobid(outlines)

    def _set_jobid(self, outlines):
        """Try to figure out the job ID from the SGE qsub output"""
        if len(outlines) > 0:
            m = re.compile(r"\d+").search(outlines[0])
            if m:
                self._jobid = int(m.group())
                for (num, worker) in enumerate(self):
                    worker._jobid = "%d.%d" % (self._jobid, num+1)
