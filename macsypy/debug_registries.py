__author__ = 'hiu-paus'

_prefix_data = '/home/hiu-paus/Projects/macsyfinder/data/'

import os


class Config(object):

    def __init__(self):
        self.profile_suffix = 'hmm'




class SystemDef(object):
    """
    Handle where are store systems. Systems are group in families and sub families.
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
        self._models = {}
        models_dir = os.path.join(self.path, 'models')
        for model in os.listdir(models_dir):
            model_path = os.path.join(models_dir, model)
            new_model = self._scan_models(model_path=model_path)
            self._models[new_model.name] = new_model


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


    def get_model(self, name):
        """
        :param name: the name of the profile to retrieve (without extension)
        :type name: string
        :returns: the absolute path of the xml definition
        :rtype: string
        :raise: KeyError if name does not match with any model definition
        """
        name_path = name.split('/')
        models = self._models
        for level in name_path:
            if level in self._models:
                model = models[level]
                models = model
            else:
                raise KeyError("{} does not match with any models".format(level))
        return model.path


    def get_profile(self, name):
        """
        :param name: the name of the profile to retrieve (without extension)
        :type name: string
        :returns: the absolute path of the hmm profile
        :rtype: string
        :raise: KeyError if name does not match with any profiles
        """
        return self._profiles[name]


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


class SystemsDefRegistry(object):
    """
    scan canonical directories to register the different systems availables
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




cfg = Config()

reg = SystemsDefRegistry(cfg)
print reg._register

cas = reg._register['CRISPR-Cas']
