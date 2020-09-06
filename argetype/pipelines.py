"""argetype pipelines module

I have put a lot of philosophical thought into what is the difference
between a pipeline and a task. In the implementation of this package
they both implement the RunInterface: they are supposed to run something
and reach a state of completion at some point.

A Task (TaskBase), when it explicitly defines all its dependencies, can 
be considered a pipeline indirectly, and it is not necessary to wrap it 
in a Pipeline class (PbBase).

A Pipeline class is useful when some glue code is necesarry between
different tasks:
    - output and input of subsequent Tasks must be mapped
    - a Task must be applied to a range of samples

Todo:
    - extracting task settings

Example:
    >>> from argetype.pipelines import PbBase, step
    ... from argetype.tasks import TaskBase
    ... 
    ... class TaskA(TaskBase):
    ...     pass
    ... 
    ... class TaskB(TaskBase):
    ...     dependence: TaskBase
    ... 
    ... class TaskC(TaskBase):
    ...     taskB: TaskB
    ... 
    ... class Pipeline(PbBase):
    ...     dependence: TaskA = STEP()
    ...     taskB: TaskB
    ...     end:   TaskC = STEP(input=('dependence', 'taskB'))
    ... 
    ... pipeline = Pipeline()
    ... pipeline.run()

"""
from .tasks import RunInterface
from itertools import count

class StepBase(object):
    lastNumber = count(0)
    
    def __init__(self, number = None, input = {}, output = {}):
        self.number = number if number else next(self.lastNumber)
        self.input = input
        self.output = output

    def __call__(self, pipeline, taskName):
        # TODO check previous ran tasks, and see if they need to be added
        task = pipeline._tasks[taskName](parse=False)
        if not task.completed():
            task(**{k:pipeline._tasks[v] for k,v in self.input.items()})
            task.run()
            pipeline._tasks[taskName] = task

STEP = StepBase

class PbBase(RunInterface):
    def __init__(self):
        import typing
        import inspect
        self._tasks = typing.get_type_hints(self)
        self.steps = list(self._tasks)
        self._step_management = {
            k:v for k,v in vars(type(self)).items()
            if k in self._tasks
        }
        for step in self.steps:
            if step not in self._step_management:
                self._step_management[step] = STEP(number=-1)
        self._ran = False

    def run(self):
        self._run_started = True
        for step in self.steps:
            self._step_management[step](
                pipeline = self,
                taskName = step
            )
        self._ran = True

    def completed(self):
        return self._ran
