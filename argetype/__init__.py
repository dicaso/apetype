"""argetype module

Provides base class `Settings` that can be used
for making a derivative class that contains typed
settings, for a module or package.
"""

import typing
from collections import OrderedDict

class ConfigBase(object):
    """Takes in a list of settings, which will be exposed
    as CLI arguments. Each settings tuple should have the
    following format:
    ('--name', keyword dict for the parser.add_argument function)

    The recommended way to build a SettingsBase object, is to
    inherit from it and define the `setup` method
    (see SettingsBase.setup docstring)

    Args:
      parse (bool): if True, already parse arguments
    """

    def __init__(self, parse=True):
        self.groups = False
        self.setup()
        self.make_parser()
        if parse:
            self.parse_args()

    def __getitem__(self, key):
        return self.settings.__getattribute__(key)

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            return self.settings.__getattribute__(name)

    def setup(self):
        """Can be overwritten by inheriting classes.
        Allows defining parameters with type hints.
        Overwritten setup methods need to call `super().setup()`
        at the end.

        Example:
        >>> class Settings(SettingsBase):
        ...     def setup(_, a: int = 5, b: float = .1, c: str = 'a'):
        ...          super().setup()
        ... settings = Settings()
        """
        import inspect
        cls = type(self)
        cls_vars = vars(cls)
        # getting hints from self does not always work
        annotations = typing.get_type_hints(cls)
        annotation_groups = [attr for attr in cls_vars if isinstance(cls_vars[attr], type)]
        if annotations or annotation_groups:
            # get_attr does not work so requesting vars
            self._settings = [
                (f'--{p}', {
                    'default': cls_vars[p],  # does not yet work for positionals
                    'type': annotations[p]
                }
                )
                for p in annotations
            ]
            if annotation_groups:
                self.groups = True
                self._settings = OrderedDict(('default', self._settings)) if self._settings else OrderedDict()
                for annot_grp in annotation_groups:
                    grp_annotations = typing.get_type_hints(cls_vars[annot_grp])
                    grp_cls_vars = vars(cls_vars[annot_grp])
                    self._settings[annot_grp] = [
                        (f'--{p}', {
                        'default': grp_cls_vars[p],  # does not yet work for positionals
                        'type': grp_annotations[p]
                        }
                        )
                        for p in grp_annotations
                    ]
        else:
            sig = inspect.signature(self.setup)
            # source = inspect.getsource(self.setup)  # to extract help comments
            # print(sig)
            self._settings = [
                (f'--{p}', {
                    'default': sig.parameters[p].default,
                    'type': sig.parameters[p].annotation
                }
                )
                for p in sig.parameters
            ]

    def make_parser(self, **kwargs):
        import argparse
        self.parser = argparse.ArgumentParser(**kwargs)
        for grp in self._settings:
            if self.groups:
                parser = self.parser.add_argument_group(grp)
            else:
                parser = self.parser
            for setting in (self._settings[grp] if self.groups else self._settings):
                parser.add_argument(setting[0], **setting[1])
            if not self.groups:
                break  # if no groups need to break

    def parse_args(self):
        self.settings = self.parser.parse_args()
        return self.settings
