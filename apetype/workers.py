"""Module for handling multiprocessing
Workers take a Task and check which subtasks
can be executed

If subtasks generate figures, matplotlib needs to
use a compatible backend. Tested backends: 'pdf'
The main process should therefore execute:

>>> import matplotlib
... matplotlib.use('pdf')

Todo:
    - connect to workers on other servers
    - distribute TaskBase depedencies over the workers

Example:
    >>> import apetype as at
    ... import apetype.workers
    ... import time, numpy as np
    ... 
    ... class Task(at.TaskBase):
    ...     a: int = 5
    ... 
    ...     def ar(_, a) -> np.ndarray:
    ...         return np.array([a*i for i in range(10)])
    ... 
    ...     def aplus10(_, ar: at.tasks.InjectItems) -> list:
    ...         print('y', ar)
    ...         time.sleep(2)
    ...         return ar+10
    ... 
    ...     def amaal10(_, ar, aplus10) -> np.ndarray:
    ...         return ar*aplus10
    ... 
    ... task = Task(parse=True)
    ... manager = at.workers.Manager(task, max_workers_x_subtask=2)
    ... manager.start()
"""
import os
import multiprocessing as mp
from .utils import CountSemaphore
from .configs import ConfigBase
from .tasks import RunInterface, InjectItems

class Manager(ConfigBase, RunInterface):
    workers: int = os.cpu_count()
    max_workers_x_subtask: int = workers
    
    def __init__(self, task, *args, **kwargs):
        # Init ConfigBase, as Manager does not have positionals parse=True
        super().__init__(parse=True)
        # Set Manager options passed as kwargs
        if args or kwargs: self(*args, **kwargs)
        self.task = task
        # task has to be instantiated, ideally also parsed
        assert hasattr(self.task, '_output') and hasattr(self.task, '_input')
        self.manager = mp.Manager()
        self.task._input = self.manager.dict()
        self.task._output = self.manager.dict()
        self.task._output_tmp = {
            k:self.manager.dict()
            for k in self.task._output_functions.keys()
            if self.hasinjectitems(k)
        }
        self.task._output_locks = {
            k: CountSemaphore(self.max_workers_x_subtask, self.manager)
            if self.hasinjectitems(k) else mp.Lock()
            for k in self.task._output_functions.keys()
        }
        self.task._output_events = {
            k: mp.Event()
            for k in self.task._output_functions.keys()
        }
        self.pool = [
            WorkerProcess(self.task, i)
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

    def hasinjectitems(self, subtask):
        for p in self.task._output_functions[subtask].parameters.values():
            if p.annotation is InjectItems: return True
        return False

    # RunInterface
    def run(self): # Usual task options not implemented
        return self.start()

    def completed(self):
        return self.task.completed()

    @property
    def _input(self):
        return self.task._input

    @property
    def _output(self):
        return self.task._output


class WorkerProcess(mp.Process):
    def __init__(self, task, number): #, queue=None):
        self.task = task
        self.number = number # used to divide subtask workload
        #self.queue = queue
        super().__init__()

    def run(self):
        for subtask in self.task._output_functions:
            lock = self.task._output_locks[subtask]
            if lock.acquire(False if isinstance(lock, mp.synchronize.Lock) else self.number):
                # Before running subtask, check if there
                # is output for all dependencies
                for dependency in self.task._output_functions[subtask].parameters.keys():
                    if dependency in self.task._output_events:
                        depevent = self.task._output_events[dependency]
                        depevent.wait()
                        # If InjectItems, register length
                        if self.task._output_functions[
                                subtask
                        ].parameters[dependency].annotation is InjectItems:
                            dynamic_injectitems_len = len(self.task._output[dependency])
                if isinstance(lock, mp.synchronize.Lock) and not subtask in self.task._output:
                    # Checking if not yet in output
                    # as due to lock racing conditions
                    # another process might have already
                    # generated the output
                    print(self.name, 'running', subtask)
                    self.task.run(subtask)
                    self.task._output_events[subtask].set()
                elif (
                        type(lock) is CountSemaphore and
                        not subtask in self.task._output and
                        self.number < dynamic_injectitems_len
                ):
                    print(self.name, 'running part of', subtask)
                    partial_result = self.task.run(
                            {subtask: range(self.number,dynamic_injectitems_len,lock._value)}, return_tmp=True
                    )
                    for k,v in partial_result.items():
                        # updating on the managed dict gave racing issues
                        self.task._output_tmp[subtask][k] = v
                    #TODO check result type in section beneath
                    # Select proc 0 to wait for subtask to finish (to avoid racing)
                    if self.number == 0 and lock._condition.acquire(False):
                        #import time
                        print(self.name, 'waiting for other processes in', subtask)
                        # Wait for all processes to finish their subtask part
                        while len(self.task._output_tmp[subtask]) < dynamic_injectitems_len:
                            lock._condition.wait(10)
                            #print(self.name, self.task._output_tmp[subtask])
                            #time.sleep(5)
                        self.task._output[subtask] = [
                            self.task._output_tmp[subtask][i]
                            for i in range(dynamic_injectitems_len)
                        ]
                        del self.task._output_tmp[subtask]
                        self.task._output_events[subtask].set()
                    else:
                        try: lock._condition.notify()
                        except RuntimeError: pass
                                
                lock.release()
