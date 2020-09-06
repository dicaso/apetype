"""argetype tasks module takes the ConfigBase
to build a inheritable TaskBase class around it.

Todo:
    - dependency settings inheritance and setting

Example:
    >>> from argtype.tasks import TaskBase
    ... class TaskDep(TaskBase):
    ...     a: str = '/tmp/file1'
    ... 
    ...     def generate_output(self) -> str:
    ...         return self.a    
    ... 
    ... class Task(TaskBase):    
    ...     # Task settings
    ...     a: int = 10
    ...     b: str = 'a'
    ... 
    ...     def generate_output1(self, task_dependence1: TaskDep) -> int:
    ...         print(task_dependence1.a)
    ...         return 0
    ...     
    ...     def generate_output2(self) -> str:
    ...         with self.env('sh') as env:
    ...             env.exec('which python')
    ...             return env.output
    ... 
    ...     def generate_output3(self) -> str:
    ...         with self.env('py') as env:
    ...             env.exec(f'''
    ...             for i in range({self.a}):
    ...                 print(i)
    ...             ''')
    ...             return env.output
    ... 
    ... task = Task()
    ... task.run()
    ... print(task._input, task._output)

""" 

import os
import abc
import typing
import inspect
from collections import OrderedDict
from argetype import ConfigBase

class RunInterface(abc.ABC):
    @abc.abstractmethod
    def run(self):
        pass

    @abc.abstractmethod
    def completed(self):
        pass

class TaskBase(ConfigBase, RunInterface):
    def __init__(self, parse=False):
        # Parsing at object creation is not that useful
        # Only when a task runs should it have all its settings.
        super().__init__(parse=parse)
        self.taskprep()
        self._input = {}
        self._output = {}
        
    def taskprep(self):
        cls = type(self)
        self._output_functions = OrderedDict([
            (name, typing.get_type_hints(fun))
            for name, fun in inspect.getmembers(cls, predicate=inspect.isfunction)
            if typing.get_type_hints(fun)
        ])

    def run(self, fail=True):
        for fn in self._output_functions:
            if fn in self._output: continue
            for dependency in self._output_functions[fn]:
                if dependency == 'return':
                    return_type = self._output_functions[fn][dependency]
                    continue
                if not dependency in self._input:
                    deptask = self._output_functions[fn][dependency]()
                    deptask.run()
                    self._input[dependency] = deptask
            return_value = self.__getattribute__(fn)(
                **{
                    dependency: self._input[dependency]
                    for dependency in self._output_functions[fn]
                    if dependency != 'return'
                }
            )
            assert isinstance(return_value, return_type)
            self._output[fn] = return_value

    def env(self, environment):
        return ExecEnvironment(environment)

    def completed(self):
        return bool(self._output)

class ExecEnvironment(object):
    predefined_envs = {
        'sh': ['bash', []],
        'py': ['python', []]
    }
    
    def __init__(self, command, options=[]):
        import threading
        import shutil
        self.lock = threading.Lock()
        self.command = shutil.which(
            self.predefined_envs[command][0]
            if command in self.predefined_envs else command
        )
        self.options = (options if options else
            self.predefined_envs[command][1]
            if command in self.predefined_envs else []
        )

    def __enter__(self):
        import tempfile
        self.lock.acquire()
        self.tmpfile = tempfile.NamedTemporaryFile(mode = 'w+t', delete = False)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # TODO log issues
        os.remove(self.tmpfile.name)
        self.lock.release()

    def exec(self, script, check = True, reident=True):
        import subprocess
        import re
        
        # reident
        if reident:
            identspace = re.compile(r'.*\n(\W*)\w')
            script = script.replace('\n'+identspace.match(script).groups()[0],'\n')
            
        # write tmp script file
        try:
            self.tmpfile.write(script)
        finally:
            self.tmpfile.close()
            
        # execute
        proc = subprocess.run(
            [self.command]+self.options+[self.tmpfile.name],
            text = True, capture_output = True, check = check
        )
        self.output = proc.stdout
        self.error = proc.stderr
