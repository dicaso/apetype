"""apetype tasks module takes the ConfigBase to build an inheritable TaskBase
class around it. Through using inheritance or including other TaskBase classes
as dependencies, the result is an implicit pipeline. For explicit pipelines
where the input and output of subsequent dependencies needs to be fitted, see
the `apetype.pipelines`.

In the example below, TaskBase dependencies is illustrated, and subtasks that
are executed in a subprocess with either a shell or python interpreter.
In the task subtasks, `self` is replaced with `_`. This is to make it clearer,
that these methods are usually not called directly by the end-user. However,
if so desired, `self` can be used. Furthermore, when writing tests, it can be
appropriate to call the subtask methods directly, allowing to control what is
being injected.

TODO:
  - for subtasks running in subprocesses, use the return_annotation to automate
    when python syntax can be used, or transformed to the syntax of the interpreter

Example:

    >>> import apetype as at
    ... class TaskDep(at.TaskBase):
    ...     a: str = '/tmp/file1'
    ... 
    ...     def generate_output(_, a) -> str:
    ...         return a    
    ... 
    ... class Task(at.TaskBase):    
    ...     # Task settings
    ...     a: int = 10
    ...     b: str = 'a'
    ... 
    ...     # Task dependencies
    ...     task_dependence1: TaskDep
    ... 
    ...     # Subtasks
    ...     def generate_output1(_, task_dependence1) -> int:
    ...         print(task_dependence1.a)
    ...         return 0
    ...     
    ...     def generate_output2(_) -> str:
    ...         with _.env('sh') as env:
    ...             env.exec('which python')
    ...             return env.output
    ... 
    ...     def generate_output3(_) -> str:
    ...         with _.env('py') as env:
    ...             env.exec(f'''
    ...             for i in range({_.a}):
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
import pickle
from collections import OrderedDict
from .configs import ConfigBase

# Interface for tasks and pipelines
class RunInterface(abc.ABC):
    @abc.abstractmethod
    def run(self):
        pass

    @abc.abstractmethod
    def completed(self):
        pass

# Interfaces for logic to perform on subtask return values
class ReturnTypeInterface(abc.ABC):
    def __init__(self, type):
        self.type = type

    @abc.abstractmethod
    def __call__(self, result):
        pass
    
    @abc.abstractmethod
    def preprocess(self):
        pass

    @abc.abstractmethod
    def postprocess(self):
        pass

class SKIP(ReturnTypeInterface):
    # SKIP subtasks
    # if subtask is explicitly mentioned in run list
    # reset type to type SKIP was instantiated with
    def __call__(self, result):
        pass
    def preprocess(self):
        # Returns SKIP instance type
        return self.type
    def postprocess(self):
        pass

class SKIPCACHE(ReturnTypeInterface):
    # Could also consider a factory method for both
    # Inheritance does not work because isinstance needs to be separate
    def __call__(self, result):
        pass
    def preprocess(self):
        # Returns SKIP instance type
        return self.type
    def postprocess(self):
        pass

# Interfaces for logic to perform on subtask value injections
class InjectInterface(abc.ABC):
    """Any InjectInterface (II) has to define a __call__ routine
    that takes the parameter that was annotated with the II, the
    subtask name and a flag dict. Flags can be modified, but 
    __call__ should only return the transformed parameter.
    """
    @abc.abstractmethod
    def __call__(self, parameter, subtask, flags):
        pass

class InjectCopy(InjectInterface):
    """To avoid having side effects on injected parameters,
    objects such as a pd.DataFrame can use this to inject
    a copy instead
    """
    def __call__(self, parameter, subtask, flags):
        return parameter.copy()

class InjectItems(InjectInterface):
    """The passed parameter should be of type list, tuple, or generator.
    An enumerate of the parameter is returned.
    """
    def __call__(self, parameter, subtask, flags):
        try:
            flags['injectitems'].add(subtask)
            assert flags['injectitems_len'] == len(parameter)
        except KeyError:
            flags['injectitems'] = {subtask}
            flags['injectitems_len'] = len(parameter)
        except AssertionError:
            raise Exception('All InjectItems params dynamic lists need to have same size')
        return zip([subtask]*len(parameter), parameter)
    
# Base class for tasks
class TaskBase(ConfigBase, RunInterface):
    def __init__(self, parse=False, prefix=True, run=False):
        """
        Parsing at object creation is not that useful
        Only when a task runs should it have all its settings.
        When run is True, parse will also be set to True.

        Args:
          parse (bool): parse settings
          prefix (bool): prefix dependency tasks or configs with attribute name
          run (bool|list): run the task upon creation, can also
            be a list of subtasks to run.
        """
        super().__init__(parse=parse or bool(run), prefix=prefix)
        self._taskprep()
        if not self.dependencies:
            # subtasks might still contain dependencies that will be
            # stored in self._input
            self._input = {}
        self._output = {}
        if run:
            self.run(
                subtasks=None if isinstance(run,bool) else run
            )
        
    def _taskprep(self):
        cls = type(self)

        # Small utility function to sort members in order of appearance
        def memberline(m):
            try:
                return m[1].__code__.co_firstlineno
            except AttributeError:
                try:
                    # Interactive code may need extra attribute
                    return m[1].__func__.__code__.co_firstlineno
                except:
                    return -1
    
        self._output_functions = OrderedDict([
            (name, inspect.signature(fun))
            for name, fun in sorted(
                    inspect.getmembers(cls, predicate=inspect.isfunction),
                    key=memberline
            )
            if typing.get_type_hints(fun)
        ])

    def _subtask_dependencies(self, subtask):
        """Get the dependencies for the subtask.
        If one of the dependencies is not available yet,
        returns False. A dependency is always first looked
        up into the task's input dict, i.c. external dependencies,
        secondly in the output dict, lastly the task attributes.

        Args:
          subtask (str): method str of the task
        
        Returns:
          dict | False
        """
        parameters = self._output_functions[subtask].parameters
        function_inputs = {}
        flags = {} # can be set by InjectInterface
        for dependency in parameters:
            if dependency in self._input:
                function_inputs[dependency] = self._input[dependency]
            # If annotation present for dependency it should have RunInterface or InjectInterface
            elif parameters[dependency].annotation is not inspect._empty:
                annotation = parameters[dependency].annotation
                if issubclass(annotation, RunInterface):
                    # Instantiate dependency and run
                    deptask = annotation() # TODO handling settings subtask
                    deptask.run()
                    self._input[dependency] = deptask
                    function_inputs[dependency] = deptask
                elif issubclass(annotation, InjectInterface):
                    try: function_inputs[dependency] = annotation()(
                            self._output[dependency], dependency, flags
                    )
                    except KeyError:
                        if dependency in self._output_functions: return False
                        else:
                            print('Currently InjectInterface only on task subtask output')
                            raise
                else:
                    raise Exception(
                        'Only RunInterface or InjectInterface supported for param annotation'
                    )
            # If no annotation it could be the output generated from one of the task methods
            elif dependency in self._output:
                function_inputs[dependency] = self._output[dependency]
            # Finally, it could also be an attribute of the task class
            else:
                try:
                    function_inputs[dependency] = self.__getattribute__(dependency)
                except AttributeError:
                    # check if attribute simply refers to 'self' or similar
                    if dependency not in ('_', 'self', 'task'):
                        if dependency in self._output_functions:
                            return False
                        else:
                            print(dependency, 'not defined/found')
                            raise
        # Return function_inputs for executing
        if not flags:
            return function_inputs
        else:
            if 'injectitems' in flags:
                parameterlist = zip(*[function_inputs.pop(k) for k in flags['injectitems']])
                return (flags['injectitems_len'], (
                    (itemnr,{**function_inputs, **dict(injitems)})
                    for itemnr, injitems in enumerate(parameterlist)
                ))
            
    def run(self, subtasks=None, fail=True, load_cache=False, return_tmp=False):
        # Task can declare a verbose attribute used in this run
        try:
            from .utils import termstyle as ts
            verbose = self.__getattribute__('verbose')
        except AttributeError: verbose = False

        # First check if external dependencies completed
        if self.dependencies:
            self._run_task_dependencies()
        
        # If subtasks is a str (when only 1 subtask needs to be executed
        # redefine it as a list to make compatible with the code
        if isinstance(subtasks, str): subtasks = [subtasks]
        
        # Run task subtask methods
        for fn in self._output_functions:
            return_annotation = self._output_functions[fn].return_annotation
            
            # If subtasks are specified, only run those
            if subtasks and fn not in subtasks: continue
            elif isinstance(return_annotation, SKIP):
                if subtasks and fn in subtasks:
                    return_annotation = return_annotation.preprocess()
                else: continue
            # If subtask already generated output continue
            if fn in self._output:
                if verbose: print(f'{ts.GREEN}{fn}{ts.RESET}', 'already generated output')
                continue
            elif hasattr(self, 'cache'):
                # If cache is set as attribute check if previous output can be loaded
                cachefilename = os.path.join(
                    self.cache,
                    f'{fn}.pickle'
                )
                if load_cache and os.path.exists(cachefilename) and not isinstance(
                        return_annotation, SKIPCACHE):
                    self._output[fn] = pickle.load(open(cachefilename, 'rb'))
                    if hasattr(self, '_cache_postprocess'):
                        # Inheriting TaskBase classes can define logic for handling cache
                        self._output[fn] = self._cache_postprocess(self._output[fn])
                    if verbose: print(f'{ts.GREEN}{fn}{ts.RESET}', 'cached output loaded')
                    continue
            return_type = return_annotation if not isinstance(
                return_annotation, SKIPCACHE) else return_annotation.preprocess()

            # TODO could check here if return type is correct in earlier generated output
            if verbose: print(f'{ts.BOLD}Executing{ts.RESET}', f'{ts.GREEN}{fn}{ts.RESET}', '...')
            function_inputs = self._subtask_dependencies(fn)
            # Checking if function_inputs is False (empty is ok = no deps)
            if function_inputs is False:
                raise Exception(f'Unmet dependencies for {fn}')
            if isinstance(function_inputs, dict):
                # regular situation -> one function_inputs
                return_value = self.__getattribute__(fn)(
                    **function_inputs
                )
            else: # function_inputs is a generator
                from collections.abc import Iterable
                complete_output, function_inputs = function_inputs
                if isinstance(subtasks, dict) and isinstance(subtasks[fn], Iterable):
                    return_value = {
                        fi[0]: self.__getattribute__(fn)(**fi[1])
                        for fi in function_inputs
                        if fi[0] in subtasks[fn]
                    }
                    if return_tmp: return return_value
                    try:
                        self._output_tmp[fn].update(return_value)
                        # Check if tmp result is complete
                        if len(self._output_tmp[fn]) == complete_output:
                            return_value = [
                                self._output_tmp[fn][i]
                                for i in sorted(self._output_tmp[fn])
                            ]
                            del self._output_tmp[fn]
                    except AttributeError:
                        self._output_tmp = {fn: return_value}
                        return
                    except KeyError:
                        self._output_tmp[fn] = return_value
                        return
                else:
                    return_value = [
                        self.__getattribute__(fn)(**fi[1])
                        for fi in function_inputs
                    ]
            if verbose: print(f'{ts.BOLD}done{ts.RESET}')
            if issubclass(return_type, ReturnTypeInterface):
                return_instance = return_type()
                self._output[fn] = return_instance(
                    task = self,
                    function = fn,
                    result = return_value
                )
            else:
                assert isinstance(return_value, return_type)
                self._output[fn] = return_value
            if hasattr(self, 'cache') and not isinstance(return_annotation, SKIPCACHE):
                # Write out output in cache dir if provided
                pickle.dump(self._output[fn], open(cachefilename, 'wb'))

    def env(self, environment):
        return ExecEnvironment(environment)

    def completed(self):
        return bool(self._output)

    # Section for dealing with external task dependencies
    def _run_task_dependencies(self):
        # Will need to be overloaded for handling multiprocessing of deps
        for dependency, deptask in self._input.items():
            if isinstance(deptask, TaskBase) and not deptask.completed():
                deptask.run()
        
class PrintInject(object):
    """Mixin class to add print diverting logic to a task.
    Classes that inherit this class, next to TaskBase can
    inject print in the subtasks.
    """
    def print(self, *args, **kwargs):
        """Method that can be used instead of print, to
        capture the stdout.

        When called without args or kwargs, returns the
        current buffer and resets it.
        """
        from io import StringIO
        if not args and not kwargs:
            printout = self._printout
            del self._printout
            return printout
        elif not 'file' in kwargs:
            out = StringIO()
            print(*args, file=out, **kwargs)
            print(out.getvalue(), end='')
            try: self._printout += out.getvalue()
            except AttributeError: self._printout = out.getvalue()
        else: print(*args, **kwargs)
        
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
            if identspace.match(script):
                # TODO should split over lines and then reassemble
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
