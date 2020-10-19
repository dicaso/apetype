import unittest

class test_Manager(unittest.TestCase):
    def setUp(self):
        from apetype.tasks import TaskBase, SKIP, SKIPCACHE, InjectCopy, InjectItems
        from apetype.workers import Manager
        
        class Task(TaskBase):
            a: int = 5
            
            def ar(_, a) -> list:
                return [a*i for i in range(10)]
            
            def aplus10(_, ar: InjectItems) -> list:
                print('y', ar)
                return ar+10
            
            def amaal10(_, ar: InjectItems, aplus10: InjectItems) -> list:
                return ar*aplus10
            
        task = Task(parse=True)
        self.manager = Manager(task, max_workers_x_subtask=2)

    def tearDown(self):
        del self.manager

    def test_manager(self):
        self.manager.start()
        self.assertEqual(len(self.manager.task._output), 3)
