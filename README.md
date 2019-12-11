# argetype

The argetype python package unites the builtin modules `argparse` and `typing`. Central is the ConfigBase class, which user classes can inherit from, to define there configurations. This is similar to configurations in the `luigi` package but with a much cleaner interface.

## Examples

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

## Perspective

This is just the initial setup of this project, but already having a basic working implementation. In the future different inherited `ConfigBase` classes should be mergeable to make a argparse parser with subparsers. The settings class could also form the basis for a task class, such as the luigi task class.

## Todo

    - write tests
    - parse comments after typed variables to serve as CLI help
    - subparser functionality
    - search a package for all defined settings classes and offer
      automated CLI interface
    - task class

