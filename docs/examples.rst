Examples
========

Configurations
--------------

    >>> from apetype import ConfigBase
    ... class Settings(ConfigBase):
    ...     positional_1: int
    ...     keyarg_1: float = .1
    ...     keyarg_2: str = 'string'
    ... settings = Settings()


