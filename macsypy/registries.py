#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2021  Institut Pasteur (Paris) and CNRS.           #
# See the COPYRIGHT file for details                                    #
#                                                                       #
# This file is part of MacSyFinder package.                             #
#                                                                       #
# MacSyFinder is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# MacSyFinder is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details .                         #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with MacSyFinder (COPYING).                                     #
# If not, see <https://www.gnu.org/licenses/>.                          #
#########################################################################

import os
from macsypy.error import MacsypyError


_separator = '/'


def split_def_name(fqn):
    """
    :param fqn: the fully qualified de name of a DefinitionLocation object
           the follow the schema model_name/<def_name>*/def_name
           for instance CRISPR-Cas/typing/cas
    :type fqn: string
    :return: the list of components of the def path
             ['CRISPR-Cas', 'typing', 'cas']
    :rtype: list of string
    """
    split = fqn.split(_separator)
    if split[0] == '':  # '/foo/bar'
        split = split[1:]
    if split[-1] == '':
        split = split[:-1]  # 'foo/bar/'
    return split


def join_def_path(*args):
    """
    join different elements of the definition path
    :param str args: the elements of the definition path, each elements must be a string
    :return: The return value is the concatenation of different elements of args with one
    separator
    :rtype: string
    """
    return _separator.join(args)


def scan_models_dir(models_dir, profile_suffix=".hmm", relative_path=False):
    """

    :param str models_dir: The path to the directory where are stored the models
    :param profile_suffix: the suffix of the hmm profiles
    :param relative_path: True if models_dir is relative false otherwise
    :return: the list of models in models_dir
    :rtype: [:class:`macsypy.registries.ModelLocation`, ...]
    """
    models = []
    for models_type in os.listdir(models_dir):
        model_path = os.path.join(models_dir, models_type)
        if os.path.isdir(model_path):
            new_model = ModelLocation(path=model_path,
                                      profile_suffix=profile_suffix,
                                      relative_path=relative_path)
            models.append(new_model)
    return models


class ModelRegistry:
    """
    scan canonical directories to register the different models available in global macsyfinder
    share data location (depending installation /usr/share/data/models) or can be
    overload with the location specify in the macsyfinder configuration (either in config file or command line)
    """

    def __init__(self):
        self._registry = {}


    def add(self, model_loc):
        """
        :param model_loc: the model location to ad to the registry
        :type model_loc: :class:`ModelLocation` object
        """
        self._registry[model_loc.name] = model_loc


    def models(self):
        """
        :returns: the list of models
        :rtype: list of :class:`ModelLocation` object
        """
        return sorted(list(self._registry.values()))  # level 0 like TXSS ou CRISPR_Cas


    def __getitem__(self, name):
        """
        :param name:
        :type name: string
        :returns: the model corresponding to name.
        :rtype: :class:`ModelLocation` object.
        :raise KeyError: if name does not match any ModelLocation registered.
        """
        if name in self._registry:
            return self._registry[name]
        else:
            raise KeyError(f"No such model definition: '{name}'")


    def __str__(self):
        s = ''

        def model_to_str(model, pad):
            if model.subdefinitions:
                model_s = f"{' ' * pad}/{model.name}\n"
                pad = pad + len(model.name) + 1
                for submodel in sorted(model.subdefinitions.values()):
                    model_s += model_to_str(submodel, pad)
            else:
                model_s = f"{' ' * pad}/{model.name}\n"
            return model_s

        for model in sorted(self.models()):
            s += model.name + '\n'
            pad = len(model.name) + 1
            for definition in model.get_definitions():
                s += model_to_str(definition, pad)
        return s


class ModelLocation:
    """
    Handle where are store Models. Models are organized in families and sub families.
    each family match to a ModelLocation. a ModelLocation contains the path toward the definitions
    and the paths to corresponding to the profiles.
    """

    def __init__(self, path=None, profile_dir=None, def_dir=None, profile_suffix='.hmm', relative_path=False):
        """
        :param str path: if it's an installed model, path is the absolute path to a model family.
                     otherwise path is None, and profile_dir and def_dir must be specified.
        :param str profile_dir: the absolute path to the directory which contains the hmm profiles files.
        :param str def_dir: The absolute path to the directory which contains the models definitions (xml files)
                            or submodels.
        :param str profile_suffix: the suffix of hmm files
        :param bool relative_path: True if you want to work with relative path, False to work with absolute path.
        :raise: MacsypyError if path is set and profile_dir or def_dir is set
        :raise: MacsypyError if profile_dir is set but not def_dir and vice versa
        """
        if path and any((profile_dir, def_dir)):
            raise MacsypyError("'path' and '{}' are incompatible arguments".format(
                'profile_dir' if profile_dir else 'def_dir'))
        elif not path and not all((profile_dir, def_dir)):
            raise MacsypyError("if 'profile_dir' is specified 'def_dir' must be specified_too and vice versa")
        self.path = path
        if path is not None:
            self.name = os.path.basename(path)
        else:
            self.name = os.path.basename(def_dir)
        if not profile_dir:
            profile_dir = os.path.join(path, 'profiles')
        self._profiles = self._scan_profiles(profile_dir,
                                             profile_suffix=profile_suffix,
                                             relative_path=relative_path)

        self._definitions = {}
        if not def_dir:
            def_dir = os.path.join(self.path, 'definitions')
            for definition in os.listdir(def_dir):
                definition_path = os.path.join(def_dir, definition)
                new_def = self._scan_definitions(def_path=definition_path)
                if new_def:  # _scan_definitions can return None if a dir is empty
                    new_def.fqn = f"{self.name}{_separator}{new_def.fqn}"
                    if new_def.subdefinitions:
                        for def_loc in new_def.subdefinitions.values():
                            def_loc.fqn = f"{self.name}{_separator}{def_loc.fqn}"
                    self._definitions[new_def.name] = new_def
        else:
            import glob
            for model_path in glob.glob(os.path.join(def_dir, '*.xml')):

                model_fqn = os.path.basename(os.path.splitext(model_path)[0])

                if not relative_path:
                    model_path = os.path.abspath(model_path)

                new_def = DefinitionLocation(name=model_fqn,
                                             path=model_path)
                self._definitions[new_def.name] = new_def


    def _scan_definitions(self, model_def=None, def_path=None):
        """
        Scan recursively the definitions tree on the file model and store
        them.

        :param model_def: the current model definition to add new submodel location
        :type model_def: :class:`DefinitionLocation`
        :param def_path: the absolute path to analyse
        :type def_path: string
        :returns: a definition location
        :rtype: :class:`DefinitionLocation` object
        """
        if os.path.isfile(def_path):
            new_def = None
            base, ext = os.path.splitext(def_path)
            if ext == '.xml':
                name = os.path.basename(base)
                new_def = DefinitionLocation(name=name,
                                             path=def_path)
            return new_def

        elif os.path.isdir(def_path):
            name = os.path.basename(def_path)
            new_def = DefinitionLocation(name=name,
                                         path=def_path)
            for model in os.listdir(def_path):
                subdef = self._scan_definitions(model_def=new_def, def_path=os.path.join(new_def.path, model))
                if subdef is not None:
                    new_def.add_subdefinition(subdef)
            return new_def


    def _scan_profiles(self, path, profile_suffix='.hmm', relative_path=False):
        """
        Store all hmm profiles associated to the model
        """
        all_profiles = {}
        for profile in os.listdir(path):
            profile_path = os.path.join(path, profile)
            if os.path.isfile(profile_path):
                base, ext = os.path.splitext(profile)
                if ext == profile_suffix:
                    all_profiles[base] = profile_path if relative_path else os.path.abspath(profile_path)
        return all_profiles

    def __lt__(self, other):
        return self.name < other.name

    def __gt__(self, other):
        return self.name > other.name

    def __eq__(self, other):
        return self.path, self.name, self._profiles, self._definitions == \
               other.path, other.name, other._profiles, other._definitions


    def get_definition(self, fqn):
        """
        :param fqn: the fully qualified name of the definition to retrieve.
                     it's complete path without extension.
                     for instance for a file with path like this:
                     models/TXSS/defintions/T3SS.xml
                     the name is: TXSS/T3SS
                     for
                     models/CRISPR-Cas/definitions/typing/CAS.xml:
                     the name is CRISPR-Cas/typing/CAS
        :type fqn: string.
        :returns: the definition corresponding to the given name.
        :rtype: a :class:`DefinitionLocation` object.
        :raise: valueError if fqn does not match with any model definition.
        """
        name_path = fqn.split(_separator)
        def_full_name = name_path[1:]
        defs = self._definitions
        definition = None
        for level in def_full_name:
            if level in defs:
                definition = defs[level]
                defs = definition.subdefinitions
            else:
                raise ValueError(f"{level} does not match with any definitions")
        return definition


    def get_all_definitions(self, root_def_name=None):
        """
        :name root_def_name: The name of the root definition to get sub definitions.
                        If root_def is None, return all definitions for this set of models
        :param root_def_name: string
        :return: the list of definitions or subdefinitions if root_def is specified for this model.
        :rtype: list of :class: DefinitionLocation` object
        :raise ValueError: if root_def_name does not match with any definitions
        """
        if root_def_name is None:
            all_defs = [def_loc for all_loc in self._definitions.values() for def_loc in all_loc.all()]
        else:
            root_def = self.get_definition(root_def_name)
            if root_def is not None:
                all_defs = root_def.all()
            else:
                raise ValueError(f"root_def_name {root_def_name} does not match with any definitions")
        return all_defs


    def get_definitions(self):
        """
        :return:
        """
        if self._definitions is not None:
            return sorted(list(self._definitions.values()))
        else:
            return {}


    def get_profile(self, name):
        """
        :param name: the name of the profile to retrieve (without extension).
        :type name: string.
        :returns: the absolute path of the hmm profile.
        :rtype: string.
        :raise: KeyError if name does not match with any profiles.
        """
        return self._profiles[name]


    def __str__(self):
        return self.name


class MetaDefLoc(type):

    @property
    def separator(cls):
        return cls._separator


class DefinitionLocation(dict, metaclass=MetaDefLoc):
    """
    Manage were definitions are stored. a Model is a xml definition and associated profiles.
    It has 3 attributes

    name: the fully qualified definitions name like TXSS/T3SS or CRISPR-cas/Typing/Cas
    path: the absolute path to the definitions or set of definitions
    subdefinitions: the subdefintions if it exists
    """

    _separator = '/'

    def __init__(self, name=None, subdefinitions=None, path=None):
        super().__init__(name=name, fqn=name, subdefinitions=subdefinitions, path=path)
        self.__dict__ = self  # allow to use dot notation to access to property here name or subdefinitions ...


    @classmethod
    def split_fqn(cls, fqn):
        return [f for f in fqn.split(cls.separator) if f]

    @classmethod
    def root_name(cls, fqn):
        return cls.split_fqn(fqn)[0]

    @property
    def family_name(self):
        return self.__class__.root_name(self.fqn)

    def __hash__(self):
        return hash((self.fqn, self.path))

    def add_subdefinition(self, subdefinition):
        """
        add new sub category of definitions to this definition

        :param subdefinition: the new definition to add as subdefinition.
        :type subdefinition: :class:`DefinitionLocation` object
        """
        if self.subdefinitions is None:
            self.subdefinitions = {}
        subdefinition.fqn = f"{self.name}{_separator}{subdefinition.fqn}"
        self.subdefinitions[subdefinition.name] = subdefinition


    def all(self):
        if not self.subdefinitions:
            return [self]
        else:
            all_leaf = []
            for definition in self.subdefinitions.values():
                for leaf in definition.all():
                    all_leaf.append(leaf)
            return all_leaf


    def __str__(self):
        return self.name

    def __eq__(self, other):
        return self.fqn, self.path, self.subdefinitions == other.fqn, other.path, other.subdefinitions

    def __lt__(self, other):
        return self.fqn < other.fqn

    def __gt__(self, other):
        return self.fqn > other.fqn
