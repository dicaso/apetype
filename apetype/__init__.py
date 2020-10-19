"""apetype module

`ConfigBase` is the starting point for anything APEtype. It is a base
class to use to turn any script into a library function or
command-line call. Simply derive a class from it, and add typed
attributes.  Those attributes will then be offered as parameters to
the call of the instances of this class, or if you pass `parse=True`
on instantiating, will be read from the command line.

To execute a program, 2 questions are important: "What should it do?"
and "What should it know to do that?". The latter is handled by `ConfigBase`,
the former is handled by the `TaskBase`, which contains the same configuration
options and derives from `ConfigBase`, but also contains *to do* logic.
Any method that has type annotation is considered a subtask, and when calling
the task `run` method all subtasks or a selection are executed in order.

Finally, in this quick start, one might want to test the outcome of
running a task and all its subtasks. `apetype.tests` defines the
`TestTask` class that can be inherited from together with, and in
second position, a TaskBase derived class. In the `TestTask` define
all the attributes (without typing) for the test run, together with
the method outcomes as simple attributes.

Example:
    >>> import apetype as at
    ... import apetype.tests
    ... class MySettings(at.ConfigBase):
    ...     pos1: int
    ...     pos2: str
    ...     kw1: int = 5
    ... settings = MySettings(parse=False)
    ... print(settings.kw1)
    ... 
    ... class MyTask(at.TaskBase, MySettings):
    ...     kw2: int = 3
    ...     def addition(_, pos1, kw2) -> int:
    ...         return pos1+kw2
    ... mytask = MyTask(parse={'pos1':1.5,'pos2':'1'})
    ... mytask.run()
    ... print(mytask._output['addition']
    ... 
    ... class MyTaskTest(at.tests.TestTask, MyTask):
    ...     pos1 = 1
    ...     pos2 = '1'
    ...     addition = 4
    ... at.tests.unittest.main()

"""

from .configs import ConfigBase
from .tasks import TaskBase
