"""apetype configs module defines the core class ConfigBase
that is also used as the starting point for Tasks and 
ReportTasks, and wrapped in Pipelines.

Settings are defined as typed attributes on an inheriting ConfigBase class

Example:

    >>> import apetype as at
    ... class DepSettings(at.ConfigBase):
    ...     b: int
    ...     b_b: float
    ...     b_c: int = 4
    ...     
    ... class Settings(at.ConfigBase):
    ...     depconfig: DepSettings
    ...     a: int
    ...     a_b: float
    ...     a_c: int = 4
    ... 
    ... settings = Settings(parse=[1,1,1,1])
    ... print(settings.depconfig_b, settings.a, settings.depconfig.b)
    1 1 1

"""

import typing
from collections import OrderedDict

class ConfigBase(object):
    """Takes in a list of settings, which will be exposed
    as CLI arguments. Each settings tuple should have the
    following format:
    ('--name', keyword dict for the parser.add_argument function)

    The recommended way to build a ConfigBase object, is to
    inherit from it and define the `_setup` method
    (see SettingsBase._setup docstring)

    Args:
      parse (bool|list|dict): if True, already parse arguments.
        Can also be a list that will then be passed on as args,
        all list items are converted to strings.
      prefix (bool): Whether or not to add the attribute name as  prefix
        to dependency settings. Prefix setting is propagated to downstream
        dependencies.
    """

    def __init__(self, parse=True, prefix=True):
        self.groups = False
        self._prefix = prefix
        self._help = self.get_arg_docs()
        if self.dependencies:
            self._set_input()
        self._setup()
        if self.dependencies:
            self._setup_from_deps()
        self._make_call()
        self._make_parser()
        if parse:
            if isinstance(parse, bool):
                self.parse_args(None)
            elif isinstance(parse, list):
                self.parse_args([str(p) for p in parse])
            elif isinstance(parse, dict):
                self(**parse)
            else:
                raise ValueError('parse was not bool|list|dict')

    def __getitem__(self, key):
        return self.settings.__getattribute__(key)

    def __getattr__(self, name):
        if 'settings' in self.__dict__ and name in self.settings:
            return self.settings.__getattribute__(name)
        else:
            # Trick from pandas to ensure jedi completion works
            return object.__getattribute__(self, name)

    def __repr__(self):
        if 'settings' in self.__dict__:
            return str(self.settings)
        else:
            return f'<Uninitialised {self.__class__.__name__}>'

    def __str__(self):
        return self.__repr__().replace(
            'Namespace(', 'Settings:\n\n').replace(', ', '\n')[:-1]

    def _setup(self):
        """Can be overwritten by inheriting classes.
        Allows defining parameters with type hints.
        Overwritten _setup methods need to call `super()._setup()`
        at the end.

        Example:
        >>> class Settings(SettingsBase):
        ...     def _setup(_, a: int = 5, b: float = .1, c: str = 'a'):
        ...          super()._setup()
        ... settings = Settings()
        """
        import inspect
        cls = type(self)
        cls_vars = {}
        for c in cls.mro()[::-1]:
            if issubclass(c, ConfigBase) and c is not ConfigBase:
                cls_vars.update(vars(c))
        # getting hints from self does not always work
        self._annotations = annotations = typing.get_type_hints(cls)
        self._annotation_groups = annotation_groups = [
            attr for attr in cls_vars if isinstance(cls_vars[attr], type)
        ]
        if annotations or annotation_groups:
            # get_attr does not work so requesting vars
            self._settings = [
                self.add_arg_format(p, cls_vars[p] if p in cls_vars else None,
                                    annotations[p], p not in cls_vars,
                                    self._help[p] if p in self._help else None
                )
                for p in annotations
                if not issubclass(annotations[p], ConfigBase)
                # in future could make exception for parseable ConfigBase class
                # also still need solution for list|set|... annotations
            ]
            if annotation_groups:
                self.groups = True
                self._settings = OrderedDict(
                    ('default', self._settings)
                ) if self._settings else OrderedDict()
                for annot_grp in annotation_groups:
                    grp_annotations = typing.get_type_hints(cls_vars[annot_grp])
                    grp_cls_vars = vars(cls_vars[annot_grp])
                    self._settings[annot_grp] = [
                        self.add_arg_format(
                            p, grp_cls_vars[p] if p in grp_cls_vars else None,
                            grp_annotations[p], p not in grp_cls_vars,
                            self._help[p] if p in self._help else None
                        )
                        for p in grp_annotations
                    ]
        else:
            sig = inspect.signature(self._setup)
            self._settings = [
                self.add_arg_format(
                    p, sig.parameters[p].default,
                    sig.parameters[p].annotation,
                    not(bool(sig.parameters[p].default)),
                    self._help[p] if p in self._help else None
                )
                for p in sig.parameters
            ]

    def _setup_from_deps(self):
        for dependency in self._input:
            if self._prefix:
                depsettings = [
                    (
                        f'{dependency}_{s[0]}',
                        (
                            f'--{dependency}-{s[1][0][2:]}' if s[1][0].startswith('--')
                            else f'{dependency}_{s[1][0]}',
                            s[1][1]
                        )
                    )
                    for s in self._input[dependency]._settings
                ]
            else:
                depsettings = self._input[dependency]._settings
            self._settings += depsettings
        
    @staticmethod
    def add_arg_format(name, default, typ, positional, help):
        if typ is bool:
            return (name, (f'--{name.replace("_","-")}', {
                    'default': default,
                    'action': 'store_const',
                    'const': not(default),
                    'help': help
            }
            ))
        else:
            return (name, (name if positional else f'--{name.replace("_","-")}', {
                    'default': default,
                    'type': typ,
                    'help': help
            }
            ))

    @classmethod
    def get_arg_docs(cls):
        import re
        import inspect
        commentdoc = re.compile(r'\s*(\w+)\s*:\s*\w+.*?#(.+)')
        multilinedoc = re.compile(
            r'\s*(?P<identifier>\w+)\s*:\s*\w+.*?\n\s*(?P<delimiter>"""|\'\'\')'
            r'(?P<annot>(.*\n)+?)\s*?(?P=delimiter)',
            re.MULTILINE)
        onelinedoc = re.compile(
            r'\s*(?P<identifier>\w+)\s*:\s*\w+.*?\n\s*(?P<delimiter>"|\')'
            r'(?P<annot>.+)(?P=delimiter)')
        try:
            argdocs = dict([
                commentdoc.match(line).groups()
                for line in inspect.getsourcelines(cls)[0] if commentdoc.match(line)
            ])
            # Clean if not identifier and strip
            argdocs = {arg:argdocs[arg].strip() for arg in argdocs if arg.isidentifier()}
            argdocs.update({
                match.group('identifier'):match.group('annot').replace('\n', ' ').strip()
                for match in onelinedoc.finditer(inspect.getsource(cls))
                if match.group('identifier').isidentifier()
            })
            argdocs.update({
                match.group('identifier'):match.group('annot').replace('\n', ' ').strip()
                for match in multilinedoc.finditer(inspect.getsource(cls))
                if match.group('identifier').isidentifier()
            })
            return argdocs
        except TypeError:
            # TypeError will be raised when class is defined interactively
            # Returning an empty dict in that case
            return {}

    def _set_settings(self, vardict):
        from argparse import Namespace
        if 'self' in vardict: vardict.pop('self')
        for attr in vardict:
            # TODO check type
            self.__setattr__(attr, vardict[attr])
            # print(attr,vardict[attr])
        # To have consistent type of self.settings using Namespace
        # to match self.parse_args
        self.settings = Namespace(**vardict)
        if self.dependencies: self._set_input_settings(vardict)

    def _make_call(self):
        from types import FunctionType
        variables = ', '.join([
            s[0] + ('' if s[1][1]['default'] is None else
                                     '='+repr(s[1][1]['default'])
            )
            for s in sorted(
                self._settings,
                key=lambda x: x[1][1]['default'] is None,
                reverse=True
            )
        ])
        call_code = compile(f'''def set_variables(self, {variables}):
            self._set_settings(locals())
        ''', "<string>", "exec")
        # Find code sub object:
        #for i in range(len(call_code.co_consts)):
        #    if type(call_code) is type(call_code.co_consts[i]):
        #        code = call_code.co_consts[i]
        #        break
        code_pos = 0 if call_code.co_consts[-1] is None else -4
        call_func = FunctionType(
            call_code.co_consts[code_pos], globals(), "set_variables",
            argdefs = call_code.co_consts[-1] if '=' in variables else None
        )
        setattr(self.__class__, '__call__', call_func)
        
    def _make_parser(self, **kwargs):
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
                parser.add_argument(
                    setting[1][0],
                    **setting[1][1]
                )
            if not self.groups:
                break  # if no groups need to break

    def parse_args(self, args):
        self.settings = self.parser.parse_args(args)
        # Overwrite default typing values
        for a in self._annotations:
            if a in self.settings:
                self.__setattr__(a, self.settings.__getattribute__(a))
        if self.dependencies: self._set_input_settings(vars(self.settings))
        return self.settings

    # Section for dealing with dependencies
    @classmethod
    def _get_dependencies(cls):
        """Recursive property that lists all dependencies.
        Returns None if no further dependencies.

        `_get_dependencies` works recursive, because the main ConfigBase object
        needs to know all the settings that have to be requested either through
        a command line or function call.
        """
        dependency_tree = []
        clsvars = cls.__annotations__
        for attr in clsvars:
            if (
                    isinstance(clsvars[attr], type) and
                    issubclass(clsvars[attr], ConfigBase)
            ):
                dependency = clsvars[attr]
                dependency_tree.append(
                    (attr, dependency, dependency._get_dependencies())
                )

        return dependency_tree if dependency_tree else None

    @property
    def dependencies(self):
        return self._get_dependencies()

    def _set_input(self):
        # Instatiate every dependency ConfigBase without parsing
        # this will be recursive as each ConfigBase will do the same for theirs
        self._input = {
            dep[0]:dep[1](parse=False, prefix=self._prefix)
            for dep in self.dependencies
        }
        # Make dependencies accesible as attribute to self ConfigBase instance
        for k,v in self._input.items():
            self.__setattr__(k, v)

    def _set_input_settings(self, vardict):
        for dependency in self._input:
            depsettingsset = {
                f'{dependency}_{s[0]}' if self._prefix else s[0]
                for s in self._input[dependency]._settings
            }
            self._input[dependency](
                **{
                    (k[len(dependency)+1:] if self._prefix else k):v
                    for k,v in vardict.items()
                    if k in depsettingsset
                }
            )
