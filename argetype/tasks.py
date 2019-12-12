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
    ... class Task:    
    ...     # Task settings
    ...     a: int = 0
    ...     b: str = 'a'
    ... 
    ...     def generate_output1(self, task_dependence1: TaskDep) -> int:
    ...         print(task_dependence1.a)
    ...         return 0
    ...     
    ...     def generate_output2(self) -> str:
    ...         return 'a'
    ... 
    ... task = Task()
    ... task.run()
    ... print(task._input, task._output)

""" 

import typing
import inspect
from collections import OrderedDict
from argetype import ConfigBase

class TaskBase(ConfigBase):
    def __init__(self):
        super().__init__()
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
