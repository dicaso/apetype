"""apetype tests Module
Defines TestTask that allows running a full test on a TaskBase derived class.
To run the test, create a class that inherits first from TestTask and then
from the Task class to be tested. Initiate and run.

The tests Module further defines all the tests for testing apetype itself.

Example:
>>> from apetype.tests import TestTask, TaskBase
... import numpy as np
... class Task(TaskBase):
...     a: int = 5
...     b: int = 3
... 
...     def ar(_, a, b) -> np.ndarray:
...         return np.array([a,2,3])+b
...     def aplus10(_, ar) -> np.ndarray:
...         return ar+10
...     def amaal10(_, ar, aplus10) -> np.ndarray:
...         return ar*aplus10
... class TestCase(TestTask, Task):
...     a = 1
...     b = 1
...     ar = np.array([3,3,4])
...     aplus10 = np.array([3,3,4])
...     amaal10 = np.array([3,3,4])
... testcase = TestCase()
... testcase.run_test()
"""
import unittest
from ..tasks import TaskBase

class TestTask(unittest.TestCase):
    """Provide the TaskBase class to be tested either
    by creating a new class that first inherits from
    TestTask and then from the TaskBase, or initiate
    a TestTask instance passing the TaskBase class as
    an argument

    Args:
      TaskClass (None|TaskBase): class to be tested.
    """
    def __init__(self, *args, TaskClass=None, **kwargs):
        # Initialize unittest.TestCase
        super().__init__(*args, **kwargs)
            
        if TaskClass or (len(self.__class__.mro()) > 3 and
                issubclass(self.__class__.mro()[3], TaskBase)):
            self._task = (
                TaskClass if TaskClass else self.__class__.mro()[3]
            )(parse=False)
            
            # Gather defined subtask outputs
            self._expected_output = {
                subtask: self.__getattribute__(subtask)
                for subtask in self._task._output_functions
            }
            
            # Gather parameter settings
            self._task_param_input = {
                p[0].replace('-',''):self.__getattribute__(p[0].replace('-',''))
                for p in self._task._settings
                if hasattr(self, p[0].replace('-',''))
            }

    def setUp(self):
        if hasattr(self, '_task'):
            self._task(**self._task_param_input)
            
    def test_run(self):
        if not hasattr(self, '_task'): return
        
        # Run the task
        self._task.run()

        # Assert task output
        for subtask in self._expected_output:
            self.assertEqual(
                self._expected_output[subtask],
                self._task._output[subtask]
            )


# Example Task and TestCase (that does not fail)
# serving as an integration test for apetype
class ExampleTask(TaskBase):
    a: int = 5
    b: int = 3
    
    def aplusb(_, a, b) -> int:
        return a+b
    def amaalaplusb(_, a, aplusb) -> int:
        return a*aplusb

class ExampleTestCase(TestTask, ExampleTask):
    a = 1
    b = 1
    aplusb = 2
    amaalaplusb = 2

