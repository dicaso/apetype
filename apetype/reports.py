"""apetype reports module

Defines a TaskReport class, that using the leopard Report class
automatically builds a report as the Task is running.

Task subtasks can inject 'print' which then ensures that any
print statements are included in the report

Example:
    >>> from apetype.reports import TaskReport, TaskSection
    ... class Report(TaskReport):
    ...     title = 'An experimental report'
    ...     outfile = '/tmp/exreport'
    ...
    ...     def section_1(_) -> TaskSection:
    ...         import pandas as pd
    ...         return [
    ...           ('tab1', pd.DataFrame({'a': [1, 2], 'b': ['c', 'd']})),
    ...           ('tab2', pd.DataFrame({'c': [1, 2], 'd': ['c', 'd']}))
    ...         ]
    ...
    ...     def section_2_with_figure(_) -> TaskSection:
    ...         import matplotlib.pyplot as plt
    ...         print('Where will this end up?')
    ...         fig, ax = plt.subplots()
    ...         ax.scatter(range(5),range(5))
    ...         return {
    ...           'figures':{'fig 1':fig},
    ...           'clearpage':True
    ...         }

"""

from .tasks import TaskBase, PrintInject, ReturnTypeInterface

class TaskReport(TaskBase, PrintInject):
    def run(self, *args, show=True, **kwargs):
        # leopard needs to have been installed
        # if not `pip install leopard`
        # args and kwargs passed to super run
        import leopard as lp
        import inspect
        self.report = lp.Report(
            title = self.title,
            author = self.author if hasattr(self, 'author') else None,
            outfile = self.outfile,
            intro = inspect.getdoc(self) or ''
        )
        # // for newlines
        self.report.settings['doubleslashnewline'] = True
        # call super run to generate output
        super().run(*args, **kwargs)
        # output the report
        self.report.outputPDF(show=show,geometry_options=True)

    def _cache_postprocess(self, previous_output):
        """_cache_postprocess defines logic to execute
        after loading cached results.

        For the ReportTask this is just adding
        cached sections to the current report.
        """
        import leopard as lp
        if isinstance(previous_output, lp.Section):
            self.report.sections.append(previous_output)
        return previous_output
        
class TaskSection(dict, ReturnTypeInterface):
    def __init__(self):
        ReturnTypeInterface.__init__(self, dict)

    def __call__(self, task, function, result):
        import inspect
        from collections import OrderedDict
        from matplotlib.figure import Figure
        import pandas as pd
        if isinstance(result, list):
            # Optionally result can be a list of tuples
            # figures and tables will then be extracted in order
            figures = OrderedDict([r for r in result if isinstance(r[1],Figure)])
            tables = OrderedDict(
                [r for r in result
                 if isinstance(r[1],pd.DataFrame)
                 or isinstance(r[1],pd.Series)]
            )
            result = dict([r for r in result if r[0] not in figures.keys() or tables.keys()])
            if figures: result['figures'] = figures
            if tables: result['tables'] = tables
        if 'title' not in result:
            result['title'] = function.replace('_', ' ').capitalize()
        if 'text' not in result and inspect.getdoc(task.__getattribute__(function)):
            result['text'] = inspect.getdoc(task.__getattribute__(function))
        if 'clearpage' not in result and hasattr(task, 'clearpage'):
            result['clearpage'] = task.clearpage
        if hasattr(task, '_printout'):
            if 'code' in result:
                result['code'] += task.print()
            else: result['code'] = task.print()
        section = self.postprocess(result)
        task.report.subs.append(section)
        return section
        
    def preprocess(self):
        pass

    def postprocess(self, result):
        import leopard as lp
        return lp.Section(**result)
