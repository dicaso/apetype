"""apetype reports module

Defines a TaskReport class, that using the leopard Report class
automatically builds a report as the Task is running.

Task subtasks can inject 'print' which then ensures that any
print statements are included in the report
"""

from .tasks import TaskBase

class TaskReport(TaskBase):
    def run(self, *args, show=True, **kwargs):
        # leopard needs to have been installed
        # if not `pip install leopard`
        # args and kwargs passed to super run
        import leopard as lp
        self.report = lp.Report(
            title = self.title,
            outfile = self.outfile
        )
        super().run(*args, **kwargs)
        report.outputPDF(show=show)

    @property
    def print(self):
        return self.report.print
