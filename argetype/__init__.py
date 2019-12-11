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
      parse (bool | list): if True, already parse arguments.
        Can also be a list that will then be passed on as args.
    """

    def __init__(self, parse=True):
        self.groups = False
        self.setup()
        self.make_parser()
        if parse:
            self.parse_args(None if isinstance(parse, bool) else parse)

    def __getitem__(self, key):
        return self.settings.__getattribute__(key)

    def __getattr__(self, name):
        if 'settings' in self.__dict__ and name in self.settings:
            return self.settings.__getattribute__(name)
        else:
            return self.__dict__[name]

    def __repr__(self):
        if 'settings' in self.__dict__:
            return str(self.settings)
        else:
            return f'<Uninitialised {self.__class__.__name__}>'

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
        self._annotations = annotations = typing.get_type_hints(cls)
        self._annotation_groups = annotation_groups = [
            attr for attr in cls_vars if isinstance(cls_vars[attr], type)
        ]
        if annotations or annotation_groups:
            # get_attr does not work so requesting vars
            self._settings = [
                self.add_arg_format(p, cls_vars[p] if p in cls_vars else None,
                                    annotations[p], p not in cls_vars
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
                        self.add_arg_format(
                            p, grp_cls_vars[p] if p in grp_cls_vars else None,
                            grp_annotations[p], p not in grp_cls_vars
                        )
                        for p in grp_annotations
                    ]
        else:
            sig = inspect.signature(self.setup)
            # source = inspect.getsource(self.setup)  # to extract help comments
            # print(sig)
            self._settings = [
                self.add_arg_format(
                    p, sig.parameters[p].default,
                    sig.parameters[p].annotation,
                    not(bool(sig.parameters[p].default))
                )
                for p in sig.parameters
            ]

    @staticmethod
    def add_arg_format(name, default, typ, positional):
        if typ is bool:
            return (f'--{name}', {
                    'default': default,
                    'action': 'store_const',
                    'const': not(default),
            }
            )
        else:
            return (name if positional else f'--{name}', {
                    'default': default,
                    'type': typ,
            }
            )
        

    def make_parser(self, **kwargs):
        import argparse
        self.parser = argparse.ArgumentParser(**kwargs)
        if self.groups: self.group_parsers = {}
        for grp in self._settings:
            if self.groups:
                parser = self.parser.add_argument_group(grp)
                self.group_parsers[grp] = parser
            else:
                parser = self.parser
            for setting in (self._settings[grp] if self.groups else self._settings):
                parser.add_argument(setting[0], **setting[1])
            if not self.groups:
                break  # if no groups need to break

    def parse_args(self, args):
        self.settings = self.parser.parse_args(args)
        # Overwrite default typing values
        for a in self._annotations:
            if a in self.settings:
                self.__setattr__(a, self.settings.__getattribute__(a))
        return self.settings
