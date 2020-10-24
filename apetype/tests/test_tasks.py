import unittest

class test_TaskBase(unittest.TestCase):
    def setUp(self):
        import apetype as at
        class DepDepTask(at.TaskBase):
            c: int
            c_b: float = 3

            def sum(_, c, c_b) -> float:
                return .5+c+c_b
            
        class DepTask(at.TaskBase):
            depdep: DepDepTask
            b: int
            b_b: float
            b_c: int = 4

            def product(_, b, b_c) -> int:
                return b*b_c
            
        class Task(at.TaskBase):
            deptask: DepTask
            a: int
            a_b: float
            a_c: int = 4

            def product(_, a, deptask) -> int:
                return a*deptask._output['product']

        self.Task = Task

    def tearDown(self):
        del self.Task

    def test_task_with_deps_wo_prefix(self):
        task = self.Task(parse=False, prefix=False)
        task.parse_args([0,0,0,0,0])
        task.run()
        self.assertEqual(task.deptask.b, 0)
        self.assertEqual(task.b, 0)
        self.assertEqual(task._output['product'], 0)
        self.assertEqual(task.deptask._output['product'], 0)
        
    def test_task_with_deps_with_prefix(self):
        task = self.Task(parse=False, prefix=True)
        task.parse_args([0,0,0,0,0])
        task.run()
        self.assertEqual(task.deptask.b, 0)
        self.assertEqual(task.deptask_b, 0)
        self.assertEqual(task._output['product'], 0)
        self.assertEqual(task.deptask._output['product'], 0)
