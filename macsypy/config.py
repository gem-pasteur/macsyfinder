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
from time import strftime
import logging
from configparser import ConfigParser, ParsingError, NoSectionError

from macsypy import __MACSY_CONF__, __MACSY_DATA__
from  macsypy.model_conf_parser import ModelConfParser
_log = logging.getLogger(__name__)


class MacsyDefaults(dict):
    """
    Handle all default values for macsyfinder.
    the default values must be defined here, **NOT** in argument parser nor in config
    the argument parser or config must use a MacsyDefaults object
    """

    def __init__(self, **kwargs):
        """
        :param kwargs: allow to overwrite a default value.
                       It mainly used in unit tests

        To define a new default value just add an attribute with the default value
        """
        super().__init__()
        self.__dict__ = self
        if __MACSY_DATA__ == '$' + 'MACSYDATA':
            prefix_data = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data'))
        else:
            prefix_data = os.path.join(__MACSY_DATA__, 'data')
        self.cfg_file = kwargs.get('cfg_file', None)
        self.coverage_profile = kwargs.get('coverage_profile', 0.5)
        self.e_value_search = kwargs.get('e_value_search', 0.1)
        self.no_cut_ga = kwargs.get('no_cut_ga', False)
        self.db_type = kwargs.get('db_type', None)
        self.hmmer = kwargs.get('hmmer', 'hmmsearch')
        self.i_evalue_sel = kwargs.get('i_evalue_sel', 0.001)
        self.idx = kwargs.get('idx', False)
        self.inter_gene_max_space = kwargs.get('inter_gene_max_space', None)
        self.log_level = kwargs.get('log_level', logging.INFO)
        self.log_file = kwargs.get('log_file', 'macsyfinder.log')
        self.max_nb_genes = kwargs.get('max_nb_genes', None)
        self.min_genes_required = kwargs.get('min_genes_required', None)
        self.min_mandatory_genes_required = kwargs.get('min_mandatory_genes_required', None)
        self.models = kwargs.get('models', [])
        self.models_dir = kwargs.get('models_dir', os.path.join(prefix_data, 'models'))
        self.multi_loci = kwargs.get('multi_loci', set())
        self.mute = kwargs.get('mute', False)
        self.out_dir = kwargs.get('out_dir', None)
        self.previous_run = kwargs.get('previous_run', None)
        self.profile_suffix = kwargs.get('profile_suffix', '.hmm')
        self.quiet = kwargs.get('quiet', 0)
        self.relative_path = kwargs.get('relative_path', False)
        self.replicon_topology = kwargs.get('replicon_topology', 'circular')
        self.res_extract_suffix = kwargs.get('res_extract_suffix', '.res_hmm_extract')
        self.res_search_dir = kwargs.get('res_search_dir', os.getcwd())
        self.res_search_suffix = kwargs.get('res_search_suffix', '.search_hmm.out')
        self.sequence_db = kwargs.get('sequence_db', None)
        self.topology_file = kwargs.get('topology_file', None)
        self.verbosity = kwargs.get('verbosity', 0)
        self.worker = kwargs.get('worker', 1)
        self.mandatory_weight = kwargs.get('mandatory_weight', 1.0)
        self.accessory_weight = kwargs.get('accessory_weight', .5)
        self.neutral_weight = kwargs.get('neutral_weight', 0.0)
        self.exchangeable_weight = kwargs.get('exchangeable_weight', .8)
        self.loner_multi_system_weight = kwargs.get('loner_multi_system_weight', .7)
        self.itself_weight = kwargs.get('itself_weight', 1.0)
        self.redundancy_penalty = kwargs.get('redundancy_penalty', 1.5)


class Config:
    """
    Handle configuration values for macsyfinder.
    This values come from default and ar superseded by the configuration files, then the command line settings.
    """

    cfg_opts = [('base', ('db_type', 'idx', 'replicon_topology', 'sequence_db', 'topology_file')),
                ('models_opt', ('inter_gene_max_space', 'max_nb_genes', 'min_mandatory_genes_required',
                                'min_genes_required', 'multi_loci')),
                ('models', tuple()),
                ('hmmer', ('coverage_profile', 'e_value_search', 'no_cut_ga', 'i_evalue_sel', 'hmmer')),
                ('score_opt', ('mandatory_weight', 'accessory_weight', 'neutral_weight', 'exchangeable_weight',
                               'itself_weight', 'redundancy_penalty', 'loner_multi_system_weight')),
                ('directories', ('models_dir', 'out_dir', 'profile_suffix', 'res_search_dir',
                                 'res_search_suffix', 'res_extract_suffix')),
                ('general', ('cfg_file', 'log_file', 'log_level', 'previous_run', 'relative_path',
                             'verbosity', 'quiet', 'mute', 'worker')),
                ]

    model_opts = ('itself', 'exchangeable', 'mandatory', 'accessory', 'neutral',
                  'loner_multi_system', 'redundancy_penalty',
                  'e_value_search', 'e_value_sel', 'coverage_profile', 'cut_ga')


    def __init__(self, defaults, parsed_args):
        """
        Store macsyfinder configuration options and propose an interface to access to them.

        The config object is populated in several steps, the rules of precedence are

        system wide conf < user home conf < model conf < (project conf | previous run) < command line

        system wide conf = etc/macsyfinder/macsyfinder.conf
        user home conf = ~/.macsyfinder/macsyfinder.conf
        model conf = model_conf.xml at the root of the model package
        project conf = macsyfinder.conf  where the analysis is run
        previous run = macsyfinder.conf in previous run results dir
        command line = the options set on the command line

        :param defaults:
        :type defaults: a :class:`MacsyDefaults` object
        :param parsed_args: the command line arguments parsed
        :type parsed_args: a :class:`argspace.Namescape` object
        """
        self.cfg_name = "macsyfinder.conf"
        self._defaults = defaults

        if __MACSY_DATA__ == '$' + 'MACSYDATA':
            self._prefix_data = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data'))
        else:
            self._prefix_data = os.path.join(__MACSY_DATA__, 'data')

        if __MACSY_CONF__ == '$' + 'MACSYCONF':
            self._conf_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'etc'))
        else:
            self._conf_dir = __MACSY_CONF__

        self._options = {}
        self._tmp_opts = {}

        self._set_default_config()

        system_wide_config_file = os.path.join(self._conf_dir, self.cfg_name)
        if os.path.exists(system_wide_config_file):
            self._set_system_wide_config(system_wide_config_file)

        user_wide_config_file = os.path.join(os.path.expanduser('~'), '.macsyfinder', self.cfg_name)
        if os.path.exists(user_wide_config_file):
            self._set_user_wide_config(user_wide_config_file)

        project_config_file = os.path.join(os.getcwd(), 'macsyfinder.conf')
        if os.path.exists(project_config_file):
            self._set_project_config_file(project_config_file)

        if hasattr(parsed_args, 'cfg_file') and parsed_args.cfg_file:
            user_config_file = parsed_args.cfg_file
            if os.path.exists(user_config_file):
                self._set_user_config_file(user_config_file)
            else:
                raise ValueError(f"Config file {user_config_file} not found.")

        previous_run = False
        if hasattr(parsed_args, 'previous_run') and parsed_args.previous_run:
            prev_config = os.path.normpath(os.path.join(parsed_args.previous_run, self.cfg_name))
            previous_run = True
            if not os.path.exists(prev_config):
                raise ValueError(f"No config file found in dir {parsed_args.previous_run}")
            self._set_previous_run_config(prev_config)

        if previous_run:
            if hasattr(parsed_args, 'sequence_db') and parsed_args.sequence_db:
                _log.warning(f"ignore sequence_db '{parsed_args.sequence_db}' use sequence_db "
                             f"from previous_run '{parsed_args.previous_run}'.")
                parsed_args.sequence_db = None

        self._set_command_line_config(parsed_args)

        models = self.models()
        if models:
            # the config is also used for macsydata
            # in this case the models are not necessary
            model_family_root, _ = models
            models_config_file = os.path.join(self.models_dir(), model_family_root, 'model_conf.xml')
            if os.path.exists(models_config_file):
                self._set_model_config(models_config_file)

        # superseed options (potentially in model_conf)
        # by the values provided by previous-run, project conf, the users on the commandline
        self._options.update(self._tmp_opts)


    def _set_options(self, options):
        """
        set key, value in the general config

        :param options: the options to specify in general config
        :type options: a dictionary with option name as keys and values as values
        """
        for opt, val in options.items():
            if val is not None:
                met_name = f'_set_{opt}'
                if hasattr(self, met_name):
                    # config has a specific method to parse and store the value
                    # for this option
                    getattr(self, met_name)(val)
                else:
                    # config has no method defined to set this option
                    self._options[opt] = val


    def _set_default_config(self):
        """
        set the value comming from MacsyDefaults
        """
        # the special methods are not used to fill with defaults values
        self._options = {k: v for k, v in self._defaults.items()}


    def _set_system_wide_config(self, config_path):
        """
        set the options from the system wide configuration file

        :param str config_path:
        """
        system_wide_config = self._config_file_2_dict(config_path)
        self._set_options(system_wide_config)


    def _set_user_wide_config(self, config_path):
        """
        Set the options from the ~/.macsyfinder/macsyfinder.conf file

        :param str config_path: The path to the ~/.macsyfinder/macsyfinder.conf
        """
        user_wide_config = self._config_file_2_dict(config_path)
        self._set_options(user_wide_config)


    def _set_model_config(self, model_conf_path):
        """
        Set the options from the model package model_conf.xml file

        :param str model_conf_path: The path to the model_conf.xml file
        """
        mp = ModelConfParser(model_conf_path)
        model_conf = mp.parse()

        self._set_options(model_conf)


    def _set_project_config_file(self, config_path):
        """
        Set the options from the macsyfinder.conf present in the current directory

        :param str config_path: the path to the configuration file
        """
        project_config = self._config_file_2_dict(config_path)
        for opt in self.model_opts:
            if opt in project_config:
                self._tmp_opts[opt] = project_config[opt]
                del(project_config[opt])
        self._set_options(project_config)


    def _set_user_config_file(self, config_path):
        """
        Set the options specified by the user on the command line via the --cfg-file option

        :param str config_path: The path to the configuration path
        """
        user_config = self._config_file_2_dict(config_path)
        for opt in self.model_opts:
            if opt in user_config:
                self._tmp_opts[opt] = user_config[opt]
                del (user_config[opt])
        self._set_options(user_config)

    def _set_previous_run_config(self, prev_config_path):
        """
        Set the options specified by the user on the command line via --previous-run

        :param prev_config_path:
        """
        previous_conf = self._config_file_2_dict(prev_config_path)
        if 'out_dir' in previous_conf:
            # set the out_dir from the previous_run is a non sense
            del(previous_conf['out_dir'])
        for opt in self.model_opts:
            if opt in previous_conf:
                self._tmp_opts[opt] = previous_conf[opt]
                del (previous_conf[opt])
        self._set_options(previous_conf)


    def _set_command_line_config(self, parsed_args):
        """

        :param parsed_args: the argument set on the command line
        :type parsed_args: :class:`argparse.Namespace` object.
        """
        # do not iter on args special attribute
        # do not set option if the value is None (not specified on command line)
        args_dict = {k: v for k, v in vars(parsed_args).items() if not k.startswith('__') and v is not None}
        for opt in self.model_opts:
            if opt in args_dict:
                self._tmp_opts[opt] = args_dict[opt]
                del (args_dict[opt])
        self._set_options(args_dict)


    def _config_file_2_dict(self, file):
        """
        Parse a configuration file in ini format in dictionnary

        :param str file: path to the configuartion file
        :return: the parsed options
        :rtype: dict
        """
        parser = ConfigParser()
        parse_meth = {int: parser.getint,
                      float: parser.getfloat,
                      bool: parser.getboolean
                      }
        try:
            parser.read([file])
            _log.debug(f"Configuration file {file} parsed.")
        except ParsingError as err:
            raise ParsingError(f"The macsyfinder configuration file '{file}' is not well formed: {err}") from None
        opts = {}
        sections = [s for s in parser.sections()]

        for section in sections:
            for option in parser.options(section):
                opt_type = type(self._defaults.get(option, None))
                try:
                    opt_value = parse_meth.get(opt_type, parser.get)(section, option)
                except (ValueError, TypeError) as err:
                    raise ValueError(f"Invalid value in config_file for option '{option}': {err}")
                opts[option] = opt_value
        return opts


    def __getattr__(self, option_name):
        # some getter return just a value they can be transformed in property
        # but some other need extra argument so they cannot be a property, they must be methods
        # to have something generic and with the same behavior
        # that mean need to call all of them
        # for generic getter, that mean no code in config
        # I simulate a function (lambda) which can be called without argument
        if option_name in self._options:
            return lambda: self._options[option_name]
        else:
            raise AttributeError(f"config object has no attribute '{option_name}'")


    def _str_2_tuple(self, value):
        """
        transform a string with syntax {model_fqn int} in list of tuple

        :param str value: the string to parse
        :return:
        :rtype: [(model_fqn, int), ...]
        """
        try:
            it = iter(value.split())
            res = [(a, next(it)) for a in it]
            return res
        except StopIteration:
            raise ValueError(f"You must provide a list of model name and value separated by spaces: {value}")


    def save(self, path_or_buf=None):
        """
        save itself in a file in ini format.

        .. note::
            the undefined options (set to None) are omitted

        :param path_or_buf: where to serialize itself.
        :type path_or_buf: str or file like object
        """
        def serialize():
            conf_str = ''
            for section, options in self.cfg_opts:
                conf_str += f"[{section}]\n"
                if section == 'models':
                    # [(model_family, (def_name1, ...)), ... ]
                    model_family, model_names = self.models()
                    conf_str += f"models = {model_family} {' '.join(model_names)}\n"
                else:
                    for opt in options:
                        opt_value = self._options[opt]
                        if opt_value is None:
                            continue
                        elif isinstance(opt_value, dict):
                            value = ""
                            for model, v in opt_value.items():
                                value += f"{model} {v} "
                            opt_value = value
                        elif isinstance(opt_value, set):
                            opt_value = ', '.join(opt_value)
                        conf_str += f"{opt} = {opt_value}\n"
            return conf_str

        if path_or_buf is None:
            path_or_buf = os.path.join(self.out_dir(), self.cfg_name)
        if isinstance(path_or_buf, str):
            with open(path_or_buf, 'w') as cfg_file:
                print(serialize(), file=cfg_file)
        else:
            print(serialize(), file=path_or_buf)


    def _set_db_type(self, value):
        """
        set value for 'db_type' option

        :param str value: the value for db_type, allowed values are :
                          'ordered_replicon', 'gembase', 'unordered'
        :raise ValueError: if value is not allowed
        """
        auth_values = ('ordered_replicon', 'gembase', 'unordered')
        if value in auth_values:
            self._options['db_type'] = value
        else:
            raise ValueError(f"db_type as unauthorized value : '{value}'.")


    def _set_inter_gene_max_space(self, value):
        """
        set value for 'inter_gene_max_space' option

        :param str value: the string parse representing the model fully qualified name
                          and it's associated value and so on
                          the model_fqn is a string, the associated value must be cast in int
        :raise ValueError: if value is not well formed
        """
        opt = {}
        if isinstance(value, str):
            try:
                value = self._str_2_tuple(value)
            except ValueError as err:
                raise ValueError(f"Invalid syntax for 'inter_gene_max_space': {err}.")
        for model_fqn, quorum in value:
            try:
                opt[model_fqn] = int(quorum)
            except ValueError:
                raise ValueError(f"The value for 'inter_gene_max_space' option for model {model_fqn} must be an integer"
                                 f", but you provided {quorum}")
        self._options['inter_gene_max_space'] = opt


    def inter_gene_max_space(self, model_fqn):
        """
        :param str model_fqn: the model fully qualifed name
        :return: the gene_max_space for the model_fqn or None if it's does not specify
        :rtype: int or None
        """
        if self._options['inter_gene_max_space']:
            return self._options['inter_gene_max_space'].get(model_fqn, None)
        else:
            return None


    def _set_max_nb_genes(self, value):
        """
        set value for 'max_nb_genes' option

        :param str value: the string parse representing the model fully qualified name
                          and it's associated value and so on
                          the model_fqn is a string, the associated value must be cast in int
        :raise ValueError: if value is not well formed
        """
        opt = {}
        if isinstance(value, str):
            try:
                value = self._str_2_tuple(value)
            except ValueError as err:
                raise ValueError(f"Invalid syntax for 'max_nb_genes': {err}.")
        for model_fqn, quorum in value:
            try:
                opt[model_fqn] = int(quorum)
            except ValueError:
                raise ValueError(f"The value for 'max_nb_genes' option for model {model_fqn} must be an integer, "
                                 f"but you provided {quorum}")
        self._options['max_nb_genes'] = opt


    def max_nb_genes(self, model_fqn):
        """
        :param str model_fqn: the model fully qualifed name
        :return: the max_nb_genes for the model_fqn or None if it's does not specify
        :rtype: int or None
        """
        if self._options['max_nb_genes']:
            return self._options['max_nb_genes'].get(model_fqn, None)
        else:
            return None


    def _set_min_genes_required(self, value):
        """
        set value for 'min_genes_required' option

        :param str value: the string parse representing the model fully qualified name
                          and it's associated value and so on
                          the model_fqn is a string, the associated value must be cast in int
        :raise ValueError: if value is not well formed
        """
        opt = {}
        if isinstance(value, str):
            try:
                value = self._str_2_tuple(value)
            except ValueError as err:
                raise ValueError(f"Invalid syntax for 'min_genes_required': {err}.")
        for model_fqn, quorum in value:
            try:
                opt[model_fqn] = int(quorum)
            except ValueError:
                raise ValueError(f"The value for 'min_genes_required' option for model {model_fqn} must be an integer, "
                                 f"but you provided {quorum}")
        self._options['min_genes_required'] = opt


    def min_genes_required(self, model_fqn):
        """
        :param str model_fqn: the model fully qualifed name
        :return: the min_genes_required for the model_fqn or None if it's does not specify
        :rtype: int or None
        """
        if self._options['min_genes_required']:
            return self._options['min_genes_required'].get(model_fqn, None)
        else:
            return None

    def _set_min_mandatory_genes_required(self, value):
        """
        set value for 'min_mandatory_genes_required' option

        :param str value: the string parse representing the model fully qualified name
                          and it's associated value and so on
                          the model_fqn is a string, the associated value must be cast in int
        :raise ValueError: if value is not well formed
        """
        opt = {}
        if isinstance(value, str):
            try:
                value = self._str_2_tuple(value)
            except ValueError as err:
                raise ValueError(f"Invalid syntax for 'min_mandatory_genes_required': {err}.")
        for model_fqn, quorum in value:
            try:
                opt[model_fqn] = int(quorum)
            except ValueError:
                raise ValueError(f"The value for 'min_mandatory_genes_required' option "
                                 f"for model {model_fqn} must be an integer, but you provided {quorum}")
        self._options['min_mandatory_genes_required'] = opt


    def min_mandatory_genes_required(self, model_fqn):
        """
        :param str model_fqn: the model fully qualifed name
        :return: the min_mandatory_genes_required for the model_fqn or None if it's does not specify
        :rtype: int or None
        """
        if self._options['min_mandatory_genes_required']:
            return self._options['min_mandatory_genes_required'].get(model_fqn, None)
        else:
            return None

    def _set_models(self, value):
        """
        :param value: The models to search as return by the command line parsing or
                      the configuration files

                      if value come from command_line
                          ['model1', 'def1', 'def2', 'def3']
                      if value come from config file
                         'set_1', 'T9SS T3SS T4SS_typeI')]
                         [(model_family, [def_name1, ...]), ... ]
        """
        if isinstance(value, str):
            # it comes from a config_file
            # value = model_family_name model1 model2
            model_family_name, *models_name = value.split()
        else:
            # it come from the command line
            # value = ['model_family_name', 'model1', 'model2']
            model_family_name = value[0]
            models_name = value[1:]
        self._options['models'] = (model_family_name, models_name)


    def out_dir(self):
        """
        :return: the path to the directory where the results are stored
        """
        out_dir = self._options['out_dir']
        if out_dir:
            return out_dir
        else:
            out_dir = os.path.join(self._options['res_search_dir'],
                                   f"macsyfinder-{strftime('%Y%m%d_%H-%M-%S')}")
            self._options['out_dir'] = out_dir
            return out_dir


    def working_dir(self):
        """
        alias to :py:meth:`config.Config.out_dir`
        """
        return self.out_dir()


    def _set_replicon_topology(self, value):
        """
        set the default replicon topology

        :param str value: 'circular' or 'linear'
        """
        auth_values = ('linear', 'circular')
        value_low = value.lower()
        new_topo = None
        for topo in auth_values:
            if topo.startswith(value_low):
                new_topo = topo
                break
        if new_topo is not None:
            self._options['replicon_topology'] = new_topo
        else:
            raise ValueError(f"replicon_topology as unauthorized value : '{value}'.")


    def _set_sequence_db(self, path):
        """
        :param str path: set the path to the sequence file (in fasta format) to analysed
        """
        if os.path.exists(path) and os.path.isfile(path):
            self._options['sequence_db'] = path
        else:
            raise ValueError(f"sequence_db '{path}' does not exists or is not a file.")


    def _set_topology_file(self, path):
        """
        test if the path exists and set it in config

        :param str path: the path to the topology file
        """
        if os.path.exists(path) and os.path.isfile(path):
            self._options['topology_file'] = path
        else:
            raise ValueError(f"topology_file '{path}' does not exists or is not a file.")


    def _set_models_dir(self, path):
        """
        :param str path: the path to the models (definitions + profiles) are stored.
        """
        if os.path.exists(path) and os.path.isdir(path):
            self._options['models_dir'] = path
        else:
            raise ValueError(f"models_dir '{path}' does not exists or is not a directory.")


    def _set_multi_loci(self, value):
        """
        :param str value: the models fqn list separated by comma of multi loc models
        """
        models_fqn = {v for v in [v.strip() for v in value.split(',')] if v}
        self._options['multi_loci'] = set(models_fqn)


    def multi_loci(self, model_fqn):
        """
        :param str model_fqn: the model fully qualified name
        :return: True if the model is multi loci, False otherwise
        :rtype: bool
        """
        return model_fqn in self._options['multi_loci']


    def hmmer_dir(self):
        """

        :return: The name of the directory containing the hmmsearch results (output, error, parsing)
        """
        return 'hmmer_results'


    def hit_weights(self):
        """

        :return: the options used in scoring systems (mandatory_weight, accessory_weight, itself_weight,
                 exchangeable_weight, loner_multi_system_weight)
        :rtype: dict
        """
        return {'mandatory': self._options['mandatory_weight'],
                'accessory': self._options['accessory_weight'],
                'neutral': self._options['neutral_weight'],
                'itself': self._options['itself_weight'],
                'exchangeable': self._options['exchangeable_weight'],
                'loner_multi_system': self._options['loner_multi_system_weight']
                }

    def log_level(self):
        """
        :return: the verbosity output level
        :rtype: int
        """
        level = self._defaults.log_level - (10 * self.verbosity()) + (10 * self.quiet())
        level = min(50, max(10, level))
        return level


class NoneConfig:

    def __getattr__(self, property):
        if property in ('multi_loci', 'min_mandatory_genes_required', 'max_nb_genes',
                        'inter_gene_max_space', 'min_genes_required'):
            return lambda x: None
        else:
            return lambda: None
