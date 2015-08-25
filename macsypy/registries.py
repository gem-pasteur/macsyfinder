# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur (Paris) and CNRS.                         #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################



import os
from macsypy.macsypy_error import MacsypyError

_prefix_data = '$PREFIXDATA'
if 'MACSY_HOME' in os.environ and os.environ['MACSY_HOME']:
    _prefix_data = os.path.join(os.environ['MACSY_HOME'], 'data')



class ModelLocation(object):
    """
    Handle where are store Models. Models are organized in families and sub families.
    each family match to a ModelLocation. a ModelLocation contains the path toward the definitions
    and the paths to corresponding to the profiles.
    """

    def __init__(self, cfg, path=None, profile_dir=None, def_dir=None):
        """
        :param cfg: the macsyfinder configuration
        :type cfg: :class:`macsypy.config.Config` object
        :param path: if it's an installed system: path is the absolute path to a system family.
                     otherwise path is None, and profile_dir and def_dir must be specified.
        :type path: string
        :param profile_dir: the absolute path to the directory which contains the hmm profiles files.
        :type profile_dir: string
        :param def_dir: the absolute path to the directory which contains the systems definitions (xml files) or submodels.
        :type def_dir: string
        :raise: MacsypyError if path is set and profile_dir or def_dir is set
        :raise: MacsypyError if profile_dir is set but not def_dir and vice versa
        """
        if path and any((profile_dir, def_dir)):
            raise MacsypyError("'path' and '{}' are incompatible arguments".format('profile_dir' if profile_dir else 'def_dir'))
        elif not path and not all((profile_dir, def_dir)):
            raise MacsypyError("if 'profile_dir' is specified 'def_dir' must be specified_too and vice versa")
        self.cfg = cfg
        self.path = path
        if path is not None:
            self.name = os.path.basename(path)
        else:
            self.name = os.path.basename(def_dir)
        if not profile_dir:
            profile_dir = os.path.join(path, 'profiles')
        self._profiles = self._scan_profiles(profile_dir)

        self._definitions = {}
        if not def_dir:
            def_dir = os.path.join(self.path, 'definitions')
            for definition in os.listdir(def_dir):
                definition_path = os.path.join(def_dir, definition)
                new_def = self._scan_definitions(def_path=definition_path)
                if new_def:  # _scan_definitions can return None if a dir is empty
                    self._definitions[new_def.name] = new_def
        else:
            import glob
            for model_path in glob.glob(os.path.join(def_dir, '*.xml')):
                new_def = DefinitionLocation(name=os.path.basename(os.path.splitext(model_path)[0]),
                                             path=os.path.abspath(model_path))
                self._definitions[new_def.name] = new_def
#

    def _scan_definitions(self, modelf_def=None, def_path=None):
        """
        Scan recursively the _definitions tree on the file system and store
        all _definitions definitions

        :param modelf_def: the current model definition to add new submodel location
        :type definition: :class:`DefinitionLocation`
        :param def_path: the absolute path to analyse
        :type def_path: string
        """
        current_dir = None
        if os.path.isfile(def_path):
            new_def = None
            base, ext = os.path.splitext(def_path)
            if ext == '.xml':
                if modelf_def is not None:
                    name = "{}/{}".format(modelf_def.name, os.path.basename(base))
                else:
                    name = os.path.basename(base)
                new_def = DefinitionLocation(name=name,
                                             path=def_path)
            return new_def

        elif os.path.isdir(def_path):
            if modelf_def is not None:
                name = "{}/{}".format(modelf_def.name, os.path.basename(def_path))
            else:
                name = os.path.basename(def_path)
            new_def = DefinitionLocation(name=name,
                                         path=def_path)
            for model in os.listdir(def_path):
                subdef = self._scan_definitions(modelf_def=new_def, def_path=os.path.join(new_def.path, model))
                if subdef is not None:
                    new_def.add_subdefinition(subdef)
            return new_def


    def _scan_profiles(self, path):
        """
        Store all hmm profiles associated to the system
        """
        all_profiles = {}
        for profile in os.listdir(path):
            profile_path = os.path.join(path, profile)
            if os.path.isfile(profile_path):
                base, ext = os.path.splitext(profile)
                if ext == self.cfg.profile_suffix:
                    all_profiles[base] = os.path.abspath(profile_path)
        return all_profiles


    def get_definition(self, name):
        """
        :param name: the name of the definition to retrieve.
                     it's complete path without extension.
                     for instance for a file with path like this:
                     systems/TXSS/defintions/T3SS.xml
                     the name is: TXSS/T3SS
                     for
                     systems/CRISPR-Cas/definitions/typing/CAS.xml:
                     the name is CRISPR-Cas/typing/CAS
        :type name: string.
        :returns: the definition corresponding to the given name.
        :rtype: a :class:`DefinitionLocation` object.
        :raise: KeyError if name does not match with any model definition.
        """
        name_path = name.split('/')
        defs = self._definitions
        for level in name_path:
            if level in defs:
                definition = defs[level]
                defs = definition.subdefinitions
            else:
                raise KeyError("{} does not match with any definitions".format(level))
        return definition


    def get_all_definitions(self, root_def_name=None):
        """
        :name root_def_name: The name of the root definition to get sub definitions.
                        If root_def is None, return all definitions for this set of models
        :param root_def_name: string
        :return: the list of definitions or subdefinitions if root_def is specified for this model.
        :rtype: list of :class: DefinitionLocation` object
        """
        if root_def_name is None:
            all_defs = [def_loc for all_loc in self._definitions.values() for def_loc in all_loc.all()]
        else:
            root_def = self.get_definition(root_def_name)
            all_defs = root_def.all()
        return all_defs


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

    def __eq__(self, other):
        return self.path == other.path and \
               self.name == other.name and \
               self._profiles == other._profiles and \
               self._definitions == other._definitions




class DefinitionLocation(dict):
    """
    Manage were models are stored. a Model is a xml definition of a system.
    It has 3 attributes
    name: the fully qualified definitions name like TXSS/T3SS or CRISPR-cas/Typing/Cas
    path: the absolute path to the definitions or set of definitions
    subdefinitions: the subdefintions if it exists
    """

    def __init__(self, name=None, subdefinitions=None, path=None):
        super(DefinitionLocation, self).__init__(name=name, subdefinitions=subdefinitions, path=path)
        self.__dict__ = self #allow to use dot notation to access to property here name or subdefinitions ...

    def add_subdefinition(self, subdefinition):
        """
        add new sub category of definitions to this definition

        :param subdefinition: the new definition to add as subdefinition.
        :type subdefinition: :class:`DefinitionLocation` object
        """
        if self.subdefinitions is None:
            self.subdefinitions = {}
        self.subdefinitions[subdefinition.name.split('/')[-1]] = subdefinition


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
        return self.name == other.name and \
               self.path == other.path and \
               self.subdefinitions == other.subdefinitions


class ModelRegistry(object):
    """
    scan canonical directories to register the different models available in global macsyfinder
    share data location (depending installation /usr/share/data/models) or can be
    overload with the location specify in the macsyfinder configuration (either in config file or command line)
    """

    def __init__(self, cfg):
        """
        :param cfg: the macsyfinder configuration
        :type cfg: :class:`macsypy.config.Config` object
        """
        self._registery = {}
        if cfg.old_data_organization():
            new_model = ModelLocation(cfg, profile_dir=cfg.profile_dir, def_dir=cfg.def_dir)
            self._registery[new_model.name] = new_model
        else:
            if cfg.models_dir:
                models_def_root = cfg.models_dir
            else:
                models_def_root = os.path.join(_prefix_data, 'macsyfinder', 'models')
            for models_type in os.listdir(models_def_root):
                model_path = os.path.join(models_def_root, models_type)
                if os.path.isdir(model_path):
                    new_model = ModelLocation(cfg, path=model_path)
                    self._registery[new_model.name] = new_model


    def models(self):
        """
        :returns: the list of models
        :rtype: list of :class:`ModelLocation` object
        """
        return self._registery.values() # level 0 like TXSS ou CRISPR_Cas


    def __str__(self):
        s = ''

        def model_to_str(model, pad):
            if model.subdefinitions:
                model_s = "{}\{}\n".format(' ' * pad, model.name)
                pad = pad + len(model.name) + 1
                for submodel in model.subdefinitions.values():
                    model_s += model_to_str(submodel, pad)
            else:
                model_s = "{}\{}\n".format(' ' * pad, model.name)
            return model_s

        for model in self.models():
            s += model.name + '\n'
            pad = len(model.name) + 1
            for definition in model.definitions:
                s += model_to_str(definition, pad)
        return s


# class ProfilesRegistry(object):
#     """
#     ProfilesRegistry register all profiles available.
#     """
#
#     def __init__(self, cfg):
#         """
#         get all profiles available in global macsyfinder share data location (depending installation /usr/share/data/profile)
#         and overload it with the location specify in the macsyfinder configuration (either in config file or command line)
#
#         :param cfg: the macsyfinder configuration
#         :type cfg: :class:`macsypy.config.Config` object
#         """
#         self._register = {}
#         global_path = os.path.join(_prefix_data, 'macsyfinder', 'systems', 'profiles')
#         self._fill_profile(global_path, cfg)
#         local_path = cfg.profile_dir
#         if local_path:
#             self._fill_profile(local_path, cfg)
#
#     def _fill_profile(self, dir_path, cfg):
#         for path in glob.glob(os.path.join(dir_path, '*' + cfg.profile_suffix)):
#             name = os.path.basename(path)
#             name = name[:-1 * len(cfg.profile_suffix)]
#             self._register[name] = path
#
#     def __getattr__(self, name):
#         return getattr(self._register, name)
#
#
# class DefinitionsRegistry(object):
#     """
#     DefinitionsRegistry register all definition systems available.
#     """
#
#
#     def __init__(self, cfg):
#         """
#         get all systems defitions available in global macsyfinder share data location ( depending installation /usr/share/data/DEF)
#         and overload it with the location specify in the macsyfinder configuration (either in config file or command line)
#
#         :param cfg: the macsyfinder configuration
#         :type cfg: :class:`macsypy.config.Config` object
#         """
#         self._register = {}
#         global_path = os.path.join(_prefix_data, 'macsyfinder', 'DEF')
#         self._fill_def(global_path)
#         local_path = cfg.def_dir
#         if local_path:
#             self._fill_def(local_path)
#
#     def _fill_def(self, dir_path):
#         for path in glob.glob(os.path.join(dir_path, '*.xml')):
#             name = os.path.basename(path)
#             name = os.path.splitext(name)[0]
#             self._register[name] = path
#
#     def __getattr__(self, name):
#         return getattr(self._register, name)
#

