"""Module for handling multiprocessing
Workers take a Task and check which subtasks
can be executed

If subtasks generate figures, matplotlib needs to
use a compatible backend. Tested backends: 'pdf'
The main process should therefore execute:

>>> import matplotlib
... matplotlib.use('pdf')

Example:
>>> from apetype.tasks import TaskBase
... from apetype.workers import Manager
... import time
... 
... class Task(TaskBase):
...     a: int = 5
... 
...     def aplus10(_, a) -> int:
...         time.sleep(60)
...         return a+10
...
...     def amaal10(_, a, aplus10) -> int:
...         return a*aplus10
... task = Task(parse=True)
... manager = Manager(task)
... manager.start()
"""
import os
import multiprocessing as mp
from .configs import ConfigBase

class Manager(ConfigBase):
    workers: int = os.cpu_count()
    
    def __init__(self, task):
        super().__init__()
        self.task = task
        # task has to be instantiated, ideally also parsed
        assert hasattr(self.task, '_output') and hasattr(self.task, '_input')
        self.manager = mp.Manager()
        self.task._input = self.manager.dict()
        self.task._output = self.manager.dict()
        self.task._output_locks = {
            k: mp.Lock()
            for k in self.task._output_functions.keys()
        }
        self.pool = [
            WorkerProcess(self.task)
            for i in range(self.workers)
        ]
        
    def start(self):
        for p in self.pool:
            p.start()
            
        # Wait for workers to finish
        for p in self.pool:
            p.join()

        # Turn input and output back to normal dicts
        self.task._input = dict(self.task._input)
        self.task._output = dict(self.task._output)
        # Close mp manager
        self.manager.shutdown()

class WorkerProcess(mp.Process):
    def __init__(self, task): #, queue=None):
        self.task = task
        #self.queue = queue
        super().__init__()

    def run(self):
        for subtask in self.task._output_functions:
            lock = self.task._output_locks[subtask]
            if lock.acquire(block=False):
                # Before running subtask, check if there
                # is output for all dependencies
                for dependency in self.task._output_functions[subtask].parameters.keys():
                    if dependency in self.task._output_locks:
                        deplock = self.task._output_locks[dependency]
                        deplock.acquire()
                        deplock.release()
                if not subtask in self.task._output:
                    # Checking if not yet in output
                    # as due to lock racing conditions
                    # another process might have already
                    # generated the output
                    print(self.name, 'running', subtask)
                    self.task.run(subtask)
                lock.release()
