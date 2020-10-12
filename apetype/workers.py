"""Module for handling multiprocessing
Workers take a Task and check which subtasks
can be executed
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
        #self.queue = mp.Queue() # using manager instead of manually through queue
        self.pool = [
            WorkerProcess(self.task)
            for i in range(self.workers)
        ]
        
    def start(self):
        for p in self.pool:
            p.start()
            
        #self.queue.put(self.task)
        #self.queue.close()
        #self.queue.join_thread()
    
        # Wait for workers to finish
        for p in self.pool:
            p.join()

class WorkerProcess(mp.Process):
    def __init__(self, task): #, queue=None):
        self.task = task
        #self.queue = queue
        super().__init__()

    def run(self):
        #task = self.queue.get() if self.queue else self.task
        self.task.run()
        self.task._output[str(self)] = 'finished'
