# APEtype
<img title="apetype logo" src="apetype_logo.svg" width="200">

The apetype python package unites the builtin modules `argparse` and
`typing`. Central is the ConfigBase class, which user classes can
inherit from, to define their configurations. This is similar to
configurations in the `luigi` package but with a much cleaner
interface. Build upon the config class is a task class, in analogy to
`luigi`.

## Documentation

https://apetype.readthedocs.io/en/latest/

## Examples
### Settings

    from apetype import ConfigBase
    class Settings(ConfigBase):
        a: int = 0 # an integer
        d: float = .1 # a float
        c: str = 'a'
    settings = Settings()
    print(settings.a, settings['a'])

This will generate a CLI with one group of arguments.

    from apetype import ConfigBase
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
    
        def generate_output(_) -> str:
            return _.a    
    
    class Task(TaskBase):    
        # Task settings
        a: int = 10
        b: str = 'a'
    
        def generate_output1(_, task_dependence1: TaskDep) -> int:
            print(task_dependence1.a)
            return 0
        
        def generate_output2(_) -> str:
            with _.env('sh') as env:
                env.exec('which python')
                return env.output
    
        def generate_output3(_) -> str:
            with _.env('py') as env:
                env.exec(f'''
                for i in range({_.a}):
                    print(i)
                ''')
                return env.output
    
    task = Task()
    task.run()
    print(task._input, task._output)


## Perspective

This is just the initial setup of this project, but already having a basic working implementation. In the future different inherited `ConfigBase` classes should be mergeable to make a argparse parser with subparsers.

## Features

    - a class that inherits from the task to test, but also from testcase mixin
        the mixin contains test data for everything that needs to be tested
      - TODO generator to create a range of testcases based on a range of attribute values
    - TODO search a package for all defined ConfigBase classes and offer
      automated CLI interface


