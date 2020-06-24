# argetype

The argetype python package unites the builtin modules `argparse` and `typing`. Central is the ConfigBase class, which user classes can inherit from, to define their configurations. This is similar to configurations in the `luigi` package but with a much cleaner interface. Build upon the config class is a task class, in analogy to `luigi`.

## Examples
### Settings

    from argetype import ConfigBase
    class Settings(ConfigBase):
        a: int = 0 # an integer
        d: float = .1 # a float
        c: str = 'a'
    settings = Settings()
    print(settings.a, settings['a'])

This will generate a CLI with one group of arguments.

    from argetype import ConfigBase
    class SettingsDeep(ConfigBase):
        class group1:
            a: int = 0
            b: float = 0.
        class group2:
            c: str = 'a'
            d: bool = 'b'
    settings = Settings()
    print(settings.a, settings['a'])

This will generate a CLI with grouped arguments, each group having the
name of the inner class.

### Tasks

    from argtype.tasks import TaskBase
    class TaskDep(TaskBase):
        a: str = '/tmp/file1'
    
        def generate_output(self) -> str:
            return self.a    
    
    class Task(TaskBase):    
        # Task settings
        a: int = 10
        b: str = 'a'
    
        def generate_output1(self, task_dependence1: TaskDep) -> int:
            print(task_dependence1.a)
            return 0
        
        def generate_output2(self) -> str:
            with self.env('sh') as env:
                env.exec('which python')
                return env.output
    
        def generate_output3(self) -> str:
            with self.env('py') as env:
                env.exec(f'''
                for i in range({self.a}):
                    print(i)
                ''')
                return env.output
    
    task = Task()
    task.run()
    print(task._input, task._output)


## Perspective

This is just the initial setup of this project, but already having a basic working implementation. In the future different inherited `ConfigBase` classes should be mergeable to make a argparse parser with subparsers.

## Todo

    - write tests (ideally also functionality to automate writing task tests)
    - parse comments after typed variables to serve as CLI help
    - subparser functionality
    - search a package for all defined settings classes and offer
      automated CLI interface
    - inherit dependency settings in task class

