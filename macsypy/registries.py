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
import glob
from .macsypy_error import MacsypyError

_prefix_data = '$PREFIXDATA'
if 'MACSY_HOME' in os.environ and os.environ['MACSY_HOME']:
    _prefix_data = os.path.join(os.environ['MACSY_HOME'], 'data')



class SystemDef(object):
    """
    Handle where are store systems. Systems are organized in families and sub families.
    each family match to a SystemDef. a System def contains the path toward the models
    and the paths to corresponding  profiles.
    """

    def __init__(self, path, cfg):
        """
        :param path: the absolute path to a system family
        :type path: string
        :param cfg: the macsyfinder configuration
        :type cfg: :class:`macsypy.config.Config` object
        """
        self.cfg = cfg
        self.path = path
        self.name = os.path.basename(path)
        self._profiles = self._scan_profiles(os.path.join(path, 'profiles'))
        self.models = {}
        models_dir = os.path.join(self.path, 'models')
        for model in os.listdir(models_dir):
            model_path = os.path.join(models_dir, model)
            new_model = self._scan_models(model_path=model_path)
            self.models[new_model.name] = new_model


    def _scan_models(self, model_def=None, model_path=None):
        """
        Scan recursively the models tree on the file system and store
        all model definition

        :param modelf_def: the current model definition to add new submodel location
        :type model_def: :class:`ModelDefLocation`
        :param model_path: the absolute path to analyse
        :type model_path: string
        """
        if os.path.isfile(model_path):
            base, ext = os.path.splitext(model_path)
            if ext == '.xml':
                new_model = ModelDefLocation(name=os.path.basename(base),
                                             path=model_path)
            return new_model
        elif os.path.isdir(model_path):
            new_model = ModelDefLocation(name=os.path.basename(model_path),
                                         path=model_path
                                         )
            for model in os.listdir(model_path):
                submodel = self._scan_models(model_def=new_model,
                                             model_path=os.path.join(new_model.path, model))
                if submodel is not None:
                    new_model.add_submodel(submodel)

            return new_model


    def _scan_profiles(self, path):
        """
        Store all hmm profiles associated to the system
        """
        all_profiles = {}
        for profile in os.listdir(path):
            if os.path.isfile(profile):
                base, ext = os.path.splitext(profile)
                if ext == self.cfg.profile_suffix:
                    all_profiles[base] = os.path.abspath(profile)
        return all_profiles




class ModelDefLocation(dict):
    """
    Manage were models are stored. a Model is a xml definition of a system.
    """

    def __init__(self, name=None, submodels=None, path=None):
        super(ModelDefLocation, self).__init__(name=name, submodels=submodels, path=path)
        self.__dict__ = self

    def add_submodel(self, submodel):
        """
        :param submodels:
        """
        if self.submodels is None:
            self.submodels = {}
        self.submodels[submodel.name] = submodel


class SystemsRegistry(object):
    """
    scan canonical directories to register the different systems available in global macsyfinder
    share data location (depending installation /usr/share/data/profile) or can be
    overload with the location specify in the macsyfinder configuration (either in config file or command line)
    """

    def __init__(self, cfg):
        """
        :param cfg: the macsyfinder configuration
        :type cfg: :class:`macsypy.config.Config` object
        """
        self._register = {}
        systems_def_root = os.path.join(_prefix_data, 'macsyfinder', 'systems')

        for systems_type in os.listdir(systems_def_root):
            system_path = os.path.join(systems_def_root, systems_type)
            if os.path.isdir(system_path):
                new_system = SystemDef(system_path, cfg)
                self._register[new_system.name] = new_system




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

