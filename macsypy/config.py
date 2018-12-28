import os
from time import strftime
import logging
from configparser import ConfigParser, ParsingError

from macsypy import __MACSY_CONF__, __MACSY_DATA__


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
        self.__dict__ = self
        if __MACSY_DATA__ == '$' + 'MACSYDATA':
            prefix_data = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data'))
        else:
            prefix_data = os.path.join(__MACSY_DATA__, 'data')
        self.cfg_file = kwargs.get('cfg_file', None)
        self.coverage_profile = kwargs.get('coverage_profile', 0.5)
        self.e_value_search = kwargs.get('e_value_search', 1.0)
        self.db_type = kwargs.get('db_type', None)
        self.hmmer = kwargs.get('hmmer', 'hmmsearch')
        self.i_evalue_sel = kwargs.get('i_evalue_sel', 0.001)
        self.idx = kwargs.get('idx', False)
        self.inter_gene_max_space = kwargs.get('inter_gene_max_space', None)
        self.log_file = kwargs.get('log_file', 'macsyfinder.log')
        self.log_level = kwargs.get('log_level', logging.WARNING)
        self.max_nb_genes = kwargs.get('max_nb_genes', None)
        self.min_genes_required = kwargs.get('min_genes_required', None)
        self.min_mandatory_genes_required = kwargs.get('min_mandatory_genes_required', None)
        self.models = kwargs.get('models', None)
        self.models_dir = kwargs.get('models_dir', os.path.join(prefix_data, 'models'))
        self.multi_loci = kwargs.get('multi_loci', False)
        self.out_dir = kwargs.get('out_dir', os.path.join(os.getcwd(),
                                                          "macsyfinder-{}".format(strftime("%Y%m%d_%H-%M-%S"))))
        self.previous_run = kwargs.get('previous_run', False)
        self.profile_suffix = kwargs.get('profile_suffix', '.hmm')
        self.relative_path = kwargs.get('relative_path', False)
        self.replicon_topology = kwargs.get('replicon_topology', 'circular')
        self.res_extract_suffix = kwargs.get('res_extract_suffix', '.res_hmm_extract')
        self.res_search_dir = kwargs.get('res_search_dir', os.getcwd())
        self.res_search_suffix = kwargs.get('res_search_suffix', '.search_hmm.out')
        self.sequence_db = kwargs.get('sequence_db', None)
        self.topology_file = kwargs.get('topology_file', None)
        self.verbosity = kwargs.get('verbosity', 0)
        self.worker = kwargs.get('worker', 1)


class Config:

    cfg_opts = [('base', ('db_type', 'idx', 'replicon_topology', 'sequence_db', 'topology_file')),
                ('models', ('inter_gene_max_space', 'max_nb_genes', 'min_mandatory_genes_required',
                            'min_genes_required', 'models', 'multi_loci')),
                ('hmmer', ('coverage_profile', 'e_value_search', 'i_evalue_sel', 'hmmer')),
                ('directories', ('models_dir', 'out_dir', 'profile_suffix', 'res_search_dir',
                                 'res_search_suffix', 'res_extract_suffix')),
                ('general', ('cfg_file', 'log_file', 'log_level', 'previous_run', 'relative_path',
                             'verbosity', 'worker'))
                ]

    def __init__(self, defaults, parsed_args):
        """
        Store macsyfinder configuration options and propose an interface to access
        to them.

        The config object is populated with the defaults then superseded with the
        value specified in configuration files and finally by the options set on the
        command line.

        :param defaults:
        :type defaults: a :class:`MacsyDefaults` object
        :param parsed_args: the command line arguments parsed
        :type parsed_args: a :class:`argspace.Namescape` object
        """
        self._cfg_name = "macsyfinder.conf"
        self._defaults = defaults

        if __MACSY_DATA__ == '$' + 'MACSYDATA':
            self._prefix_data = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data'))
        else:
            self._prefix_data = os.path.join(__MACSY_DATA__, 'data')

        if __MACSY_CONF__ == '$' + 'MACSYCONF':
            self._conf_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),
                                                           '..', 'etc'
                                                           ))
        else:
            self._conf_dir = __MACSY_CONF__

        if hasattr(parsed_args, 'previous_run') and parsed_args.previous_run:
            prev_config = os.path.normpath(os.path.join(parsed_args.previous_run,
                                                        self._cfg_name))
            if not os.path.exists(prev_config):
                raise ValueError("No config file found in dir {}".format(parsed_args.previous_run))
            config_files = [prev_config]
        elif hasattr(parsed_args, 'cfg_file') and parsed_args.cfg_file:
            config_files = [parsed_args.cfg_file]
        else:
            config_files = [os.path.join(self._conf_dir, self._cfg_name),
                            os.path.join(os.path.expanduser('~'), '.macsyfinder', self._cfg_name),
                            'macsyfinder.conf']

        config_files_values = self._config_file_2_dict(defaults, config_files)
        args_dict = {k: v for k, v in vars(parsed_args).items() if not k.startswith('__')}

        # the special methods are not used to fill with defaults values
        self._options = {k: v for k, v in defaults.items()}

        for bag_of_opts in config_files_values, args_dict:
            for opt, val in bag_of_opts.items():
                if val is not None:
                    met_name = '_set_{}'.format(opt)
                    if hasattr(self, met_name):
                        # config has a specific method to parse and store the value
                        # for this option
                        getattr(self, met_name)(val)
                    else:
                        # config has no method defined to set this option
                        self._options[opt] = val


    def __getattr__(self, option_name):
        # some getter return just a value they can be transformed in property
        # but some other need extra argument so they cannot be a property, they must be methods
        # to have something generic and with the same behavior
        # that mean need to call all of them
        # for generic getter, that mean no code in config
        # I simulate a function (lambda) which can be called without argument
        if option_name in self._options:
            return lambda : self._options[option_name]
        else:
            raise AttributeError("config object has no attribute '{}'".format(option_name))


    def _str_2_tuple(self, value):
        """
        transform a string with syntax  {model_fqn int} n in list of tuple
        :param str value: the string to parse
        :return:
        :rtype: [(model_fqn, int), ...]
        """
        try:
            it = iter(value.split())
            res = [(a, next(it)) for a in it]
            return res
        except StopIteration:
            raise ValueError("You must provide a list of model name and"
                             " value separated by spaces: {}".format(value))


    def _config_file_2_dict(self, defaults, files):
        """
        parse config files files, the last one have precedence on the previous on so on, and return a dict
        with properties, values.
        The defaults is just used to know the type of the properties and cast them. It is not used to fill
        the dict with default values.

        :param defaults: the macsyfinder defaults value
        :type defaults: a :class:`macsypy.config.MacsyDefaults` object
        :param files: the configuration files to parse
        :type files: list of string
        :return: dict
        """
        parser = ConfigParser()
        parse_meth = {int: parser.getint,
                      float: parser.getfloat,
                      bool: parser.getboolean
                      }
        try:
            used_files = parser.read(files)
        except ParsingError as err:
            raise ParsingError("A macsyfinder configuration file is not well formed: {}".format(err)) from None

        opts = {}
        for section in parser.sections():
            for option in parser.options(section):
                opt_type = type(defaults.get(option, None))
                opt_value = parse_meth.get(opt_type, parser.get)(section, option)
                opts[option] = opt_value
        return opts

    def save(self, path_or_buf=None):
        """
        save itself in a file in ini format.

        .. note::
            the undefined options (set to None) are omitted

        :param path_or_buf: where to serialize itself.
        :type path_or_buf: str or file like object
        """
        def save(fh):
            for section, options in self.cfg_opts:
                print("[{}]".format(section), file=fh)
                for opt in options:
                    opt_value = self._options[opt]
                    if opt_value is None:
                        continue
                    elif isinstance(opt_value, dict):
                        value = ""
                        for model, v in opt_value.items():
                            value += "{} {} ".format(model, v)
                        opt_value = value
                    print("{} = {}".format(opt, opt_value), file=fh)

        if path_or_buf is None:
            path_or_buf = os.path.join(self.out_dir(), self._cfg_name)
        if isinstance(path_or_buf, str):
            with open(path_or_buf, 'w') as cfg_file:
                save(cfg_file)
        else:
            save(path_or_buf)

    def _set_db_type(self, value):
        """
        set value for 'db_type' option
        :param str value: the value for db_type, allowed values are :
                          'unordered_replicon', 'ordered_replicon', 'gembase', 'unordered'
        :raise ValueError: if value is not allowed
        """
        auth_values = ('unordered_replicon', 'ordered_replicon', 'gembase', 'unordered')
        if value in auth_values:
            self._options['db_type'] = value
        else:
            raise ValueError("db_type as unauthorized value : '{}'.".format(value))


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
                raise ValueError("Invalid syntax for 'inter_gene_max_space': {}.".format(err))
        for model_fqn, quorum in value:
            try:
                opt[model_fqn] = int(quorum)
            except ValueError:
                raise ValueError("The value for 'inter_gene_max_space' option for model {} must be an integer, "
                                 "but you provided {}".format(model_fqn, quorum))
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
                raise ValueError("Invalid syntax for 'max_nb_genes': {}.".format(err))
        for model_fqn, quorum in value:
            try:
                opt[model_fqn] = int(quorum)
            except ValueError:
                raise ValueError("The value for 'max_nb_genes' option for model {} must be an integer, "
                                 "but you provided {}".format(model_fqn, quorum))
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
                raise ValueError("Invalid syntax for 'min_genes_required': {}.".format(err))
        for model_fqn, quorum in value:
            try:
                opt[model_fqn] = int(quorum)
            except ValueError:
                raise ValueError("The value for 'min_genes_required' option for model {} must be an integer, "
                                 "but you provided {}".format(model_fqn, quorum))
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
                raise ValueError("Invalid syntax for 'min_mandatory_genes_required': {}.".format(err))
        for model_fqn, quorum in value:
            try:
                opt[model_fqn] = int(quorum)
            except ValueError:
                raise ValueError("The value for 'min_mandatory_genes_required' option "
                                 "for model {} must be an integer, but you provided {}".format(model_fqn, quorum))
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
        pass


    def out_dir(self):
        out_dir = self._options['out_dir']
        return out_dir

    def working_dir(self):
        return self.out_dir()


    def _set_replicon_topology(self, value):
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
            raise ValueError("replicon_topology as unauthorized value : '{}'.".format(value))


    def _set_sequence_db(self, path):
        if os.path.exists(path) and os.path.isfile(path):
            self._options['sequence_db'] = path
        else:
            raise ValueError("sequence_db '{}' does not exists or is not a file.".format(path))


    def _set_topology_file(self, path):
        if os.path.exists(path) and os.path.isfile(path):
            self._options['topology_file'] = path
        else:
            raise ValueError("topology_file '{}' does not exists or is not a file.".format(path))


    def _set_models_dir(self, path):
        """

        :param path:
        :return:
        """
        if os.path.exists(path) and os.path.isdir(path):
            self._options['models_dir'] = path
        else:
            raise ValueError("models_dir '{}' does not exists or is not a directory.".format(path))
