##########################################################################
#  MacSyFinder - Detection of macromolecular systems in protein dataset  #
#                using systems modelling and similarity search.          #
#  Authors: Sophie Abby, Bertrand Neron                                  #
#  Copyright (c) 2014-2022  Institut Pasteur (Paris) and CNRS.           #
#  See the COPYRIGHT file for details                                    #
#                                                                        #
#  This file is part of MacSyFinder package.                             #
#                                                                        #
#  MacSyFinder is free software: you can redistribute it and/or modify   #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  MacSyFinder is distributed in the hope that it will be useful,        #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details .                         #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with MacSyFinder (COPYING).                                     #
#  If not, see <https://www.gnu.org/licenses/>.                          #
##########################################################################
"""
Entrypoint for macsyconfig command
which generate a MacSyFinder config file
"""


import itertools
import sys
import os.path
import argparse
import shutil
from dataclasses import dataclass
from configparser import ConfigParser

from colorama import Fore, Style
from colorama import init as col_init

from macsypy import __version__ as msf_vers
from macsypy.config import MacsyDefaults
from macsypy.error import MacsypyError


class ConfigParserWithComments(ConfigParser):
    """
    Extend ConfigParser to allow comment in serialization
    """

    def add_comment(self, section, option, comment,
                    comment_nb=itertools.count(1),
                    add_space_before=False,
                    add_space_after=True):
        """
        Write a comment in .ini-format (start line with #)

        :param section: the name of the sction
        :param str option: the name of the option
        :param str comment: the comment linked to this option
        :param int comment_nb: the identifier of the comment by default an integer
        :param bool add_space_before:
        :param bool add_space_after:
        """
        comment = ''.join([f"# {l}\n" for l in comment.split('\n')])
        if add_space_before:
            comment = '\n' + comment
        if add_space_after:
            comment = comment + '\n'
        self.set(section, f"{option}_{next(comment_nb)}_comment", comment)

    def write(self, file):
        """
        Write an .ini-format representation of the configuration state.

        :param file file: the file object wher to write the configuration
        """
        for section in self._sections:
            file.write(f"[{section}]\n")
            for key, value in self._sections[section].items():
                self._write_item(file, key, value)
        file.write("\n")


    def _write_item(self, file, key, value):
        if key.endswith('_comment'):
            file.write(f"{value}")
        else:
            file.write(f"{key} = {value}\n")


@dataclass(frozen=True)
class Theme:
    """Handle color combination to hylight interactive question"""

    ERROR: str = Style.BRIGHT + Fore.RED
    WARN: str = Fore.YELLOW
    SECTION: str = Fore.MAGENTA
    RESET: str = Style.RESET_ALL
    RETRY: str = Fore.YELLOW
    QUESTION: str = Fore.GREEN
    EMPHASIZE: str = Style.BRIGHT
    EXPLANATION: str = Style.RESET_ALL
    DEFAULT: str = Style.BRIGHT + Fore.GREEN

# default theme = dark bg
theme = Theme()


def _validator(cast_func, raw, default, sequence=False):
    if raw == '':
        if default is None:
            raise MacsypyError('Please enter some value')
        else:
            raw = default
    elif sequence:
        raw = [item.strip() for item in raw.split(',')]

    try:
        if isinstance(raw, type([])):
            value = [cast_func(item) for item in raw]
        else:
            value = cast_func(raw)
    except ValueError as err:
        raise MacsypyError(f'Invalid value: {err}') from err
    return value


def check_exe(raw, default, expected, sequence=False):
    """
    Check if value point to an executable

    :param str raw: the value return by the user
    :param str default: the default value for the option
    :param expected: not used here to have the same signature for all check_xxx functions
    :return: value
    :raise MacsypyError: if the value cannot be cast in right type
    """
    def exe(value):
        exe = shutil.which(value)
        if exe:
            return value
        else:
            raise ValueError(f"'{value}' NO executable found")
    return _validator(exe, raw, default, sequence=sequence)


def check_positive_int(raw, default, expected, sequence=False):
    """
    Check if value can be cast in integer >=0

    :param str raw: the value return by the user
    :param int default: the default value for the option
    :param expected: not used here to have the same signature for all check_xxx functions
    :return: value
    :raise MacsypyError: if the value cannot be cast in right type
    """
    def positive_int(value):
        casted = int(str(value))
        if casted < 0:
            raise ValueError(f"'{value}' is not >=0")
        return casted
    return _validator(positive_int, raw, default, sequence=sequence)


def check_float(raw, default, expected, sequence=False):
    """
    Check if value can be cast in float

    :param str raw: the value return by the user
    :param float default: the default value for the option
    :param expected: not used here to have the same signature for all check_xxx functions
    :return: value
    :raise MacsypyError: if the value cannot be cast in right type
    """
    return _validator(float, raw, default, sequence=sequence)


def check_str(raw, default, expected, sequence=False):
    """
    Check if value can be cast in str

    :param str raw: the value return by the user
    :param str default: the default value for the option
    :param expected: not used here to have the same signature for all check_xxx functions
    :return: value
    :raise MacsypyError: if the value cannot be cast in right type
    """
    return _validator(str, raw, default, sequence=sequence)


def check_bool(raw, default, expected, sequence=False):
    """
    Check if value can be cast in str

    :param str raw: the value return by the user
    :param str default: the default value for the option
    :param expected: not used here to have the same signature for all check_xxx functions
    :return: value
    :raise MacsypyError: if the value cannot be cast in right type
    """
    def bool_cast(raw):
        raw = str(raw).lower()
        if raw in ('false', 'no', '0'):
            casted = False
        elif raw in ('true', 'yes', '1'):
            casted = True
        else:
            raise ValueError("Authorized values ['True'/False/0/1]")
        return casted
    return _validator(bool_cast, raw, default, sequence=sequence)


def check_dir(raw, default, expected, sequence=False):
    """
    Check if value point to a directory

    :param str raw: the value return by the user
    :param str default: the default value for the option
    :param expected: not used here to have the same signature for all check_xxx functions
    :return: value
    :raise MacsypyError: if the value cannot be cast in right type
    """
    def path(value):
        if os.path.exists(value):
            if os.path.isdir(value):
                return value
            else:
                raise ValueError(f"'{value}' is not a directory.")
        else:
            raise ValueError(f"'{value}' no such file or directory.")
    return _validator(path, raw, default, sequence=sequence)


def check_file(raw, default, expected, sequence=False):
    """
    Check if value point to a file

    :param str raw: the value return by the user
    :param str default: the default value for the option
    :param expected: not used here to have the same signature for all check_xxx functions
    :return: value
    :raise MacsypyError: if the value cannot be cast in right type
    """
    def path(value):
        if value.lower() == 'none':
            return None
        if os.path.exists(value):
            if os.path.isfile(value):
                return value
            else:
                raise ValueError(f"'{value}' is not a file.")
        else:
            raise ValueError(f"'{value}' no such file or directory.")
    return _validator(path, raw, default, sequence=sequence)


def check_choice(raw, default, expected, sequence=False):
    """
    Check if value is in list of expected values

    :param str raw: the value return by the user
    :param str default: the default value for the option
    :param expected: the allowed vlaues for this option
    :return: value
    :raise MacsypyError: if the value cannot be cast in right type
    """
    def isin(value):
        if value not in expected:
            raise ValueError(f"Authorized values are {expected}.")
        if value.lower() == 'none':
            value = None
        return value
    return _validator(isin, raw, default, sequence=sequence)


def ask(question, validator, default=None, expected=None,
        explanation='',
        sequence=False,
        question_color=None,
        retry=2):
    """
    ask a question on the terminal and return the user response
    check if the user response is allowed (right type, among allowed values, ...)

    :param str question: The question to prompt to the user on the terminal
    :param validator: what validator to be used to check the user response
    :type validator: a function define in this module starting by check\_
    :param default: the default value
    :param expected: the values allowed (can be a list of value
    :param str explanation: some explanation about the option
    :param bool sequence: True if the parameter accept a sequence of value (comma separated values)
    :param question_color: the color of the question display to the user
    :type question_color: an attribute of :class:`macsypy.scripts.macsyconfig.Theme`
    :param int retry: The number of time to repeat the question if the response is rejected
    :return: the value casted in right type
    """
    if question_color is None:
        question_color = theme.QUESTION

    question_mark = f"{question_color}?{theme.RESET}"
    if default is not None:
        if isinstance(default, type([])):
            default_str = ', '.join([str(item) for item in default])
        else:
            default_str = str(default)
        default_formatted = f" [{theme.DEFAULT}{default_str}{theme.RESET}]"
    else:
        default_formatted = ''
    if expected:
        space = '' if explanation else ' '
        expected_formatted = f"{space}{theme.QUESTION}({'/'.join([str(i) for i in expected])}){theme.RESET}"
    else:
        expected_formatted = ''
    if explanation:
        formatted_explanation = f"\n{theme.EXPLANATION}{explanation}{theme.RESET}\n"
    else:
        formatted_explanation = ''
    formatted_question = f"{question_color}> {question}{theme.RESET}"
    if explanation:
        formatted_question = formatted_question + question_mark

    raw = input(f"{formatted_question}{formatted_explanation}{expected_formatted}{default_formatted}{question_mark} ")

    try:
        val = validator(raw, default, expected, sequence=sequence)
    except MacsypyError as err:
        print(err)
        if retry > 0:
            print(f"{theme.RETRY}* {err}{theme.RESET}")
            return ask(question, validator, default=default, expected=expected, retry=retry -1)
        else:
            raise RuntimeError(f'{theme.ERROR}Too many error. Exiting{theme.RESET}') from None
    return val


def set_section(sec_name, options, config, defaults, use_defaults=False):
    """
    iter over options of a section
    ask question for each option
    and set this option in the config

    :param str sec_name: the name of the section
    :param dict options: a dictionnary with the options to set up for this section
    :param config: The config to fill in.
    :type config: :class:`ConfigParserWithComments` object
    :param defaults: the macsyfinder defaults values
    :type defaults: :class:`macsypy.config.MacsyDefaults` object
    :param bool use_defaults: The user skip this section so use defaults to set in config object
    :return:
    """

    config.add_section(sec_name)
    print(f"{theme.SECTION}Configuring {sec_name} options:{theme.RESET}\n")
    for opt_name in options:
        option = options[opt_name]
        config.add_comment(sec_name, opt_name, option['question'],
                           add_space_before=True, add_space_after=False)
        if option['explanation']:
            space_after = 'expected' in option
            config.add_comment(sec_name, opt_name, option['explanation'],
                               add_space_before=False, add_space_after=space_after)
        if 'expected' in option:
            expected = f"[{'/'.join([str(i) for i in option['expected']])}]"
            config.add_comment(sec_name, opt_name, expected,
                               add_space_before=False, add_space_after=True)
        sequence = 'sequence' in option and option['sequence']

        if use_defaults:
            value = defaults[opt_name]
        else:
            value = ask(option['question'], option['validator'],
                        default=option['default'],
                        explanation=option['explanation'],
                        expected=option.get('expected', None),
                        sequence=sequence)
        if value == defaults[opt_name]:
            if isinstance(value, type([])):
                value = ', '.join([str(item) for item in value])
            config.add_comment(sec_name, opt_name, f"{opt_name} = {value}", add_space_before=False, add_space_after=True)
        else:
            if isinstance(value, type([])):
                config.set(sec_name, opt_name, ', '.join([str(item) for item in value]))
            else:
                config.set(sec_name, opt_name, str(value))
            print()

    return config


def set_path_options(config, defaults, use_defaults=False):
    """
    Options for directories section

    :param config: The config to setup
    :type config: :class:`ConfigParserWithComments` object
    :param defaults: the macsyfinder defaults values
    :type defaults: :class:`macsypy.config.MacsyDefaults` object
    :param bool use_defaults: If True do not ask any question use the defaults values
    """
    options = {'system_models_dir': {'question': "The directory where to store the models",
                                     'validator': check_dir,
                                     'default': defaults.system_models_dir,
                                     'explanation':
"""This directory will be used as default but could be overwritten on the command line.
It will be used by macsydata to install models and macsyfinder to find them.
MacSyFinder will look for models in these directories:
 - '/share/macsyfinder/models', '/usr/local/share/macsyfinder/models'
 or
 - in ${VIRTUAL_ENV}/share/macsyfinder/models
 or
 - values provided specified by macsyfinder.conf file

then in $HOME/.macsyfinder/models and in command line option --models-dir.""",
                                     'sequence': True},
               'res_search_dir': {'question': "Results research directory",
                                  'validator': check_dir,
                                  'default': ".",
                                  'explanation':
"""macsyfinder generate a directory with all results for each jobs.
this option specify where to create these directories."""},

               'res_search_suffix': {'question': "The suffix of hmmer output",
                                     'validator': check_str,
                                     'default': defaults.res_search_suffix,
                                     'explanation': ""},

               'res_extract_suffix': {'question': "The suffix of the hmmer parsed by macsyfinder",
                                      'validator': check_str,
                                      'default': defaults.res_extract_suffix,
                                      'explanation': ""},

               'profile_suffix': {'question': "The suffix of profiles",
                                  'validator': check_str,
                                  'default': defaults.profile_suffix,
                                  'explanation': "The HMM profile provides with the models"},
               }
    if not use_defaults:
        enter = ask("Do you want to enter path options section?",
                    check_choice,
                    expected=["Y", "n"],
                    default="Y",
                    explanation="where are models, default file suffix, ...",
                    question_color=theme.EMPHASIZE + theme.SECTION
                    )
        use_defaults = enter == "no"
    set_section('directories', options, config, defaults, use_defaults=use_defaults)


def set_hmmer_options(config, defaults, use_defaults=False):
    """
    Options for hmmer section

    :param config: The config to setup
    :type config: :class:`ConfigParserWithComments` object
    :param defaults: the macsyfinder defaults values
    :type defaults: :class:`macsypy.config.MacsyDefaults` object
    :param bool use_defaults: If True do not ask any question use the defaults values
    """
    options = {'hmmer': {'question': "The binary used to search the data bank with the profiles.",
                         'validator': check_exe,
                         'default': defaults.hmmer,
                         'explanation': """If hmmer is set to None, it means that 'hmmsearch' is not found on PATH.
Ensure that 'hmmsearch' will be on the PATH at runtime or specify the 'hmmsearch' path here."""},
               'cut_ga': {'question': "Use the GA score when search with hmmsearch",
                             'validator': check_bool,
                             'default': 'Yes',
                             'expected': ['Yes', 'No'],
                             'explanation':
"""By default MSF try to applied a threshold per profile by using the
hmmer -cut-ga option. This is possible only if the GA bit score is present in the profile otherwise
MSF switch to use the --e-value-search (-E in hmmsearch).
If this option is not set the --e-value-search option is used for all profiles regardless the presence of
the a GA bit score in the profiles."""},

               'e_value_search': {'question': "Maximal e-value for hits to be reported during hmmsearch search.",
                                  'validator': check_float,
                                  'default': defaults.e_value_search,
                                  'explanation':
"""By default MSF set per profile threshold for hmmsearch run (hmmsearch --cut_ga option)
for profiles containing the GA bit score threshold.
If a profile does not contains the GA bit score the --e-value-search (-E in hmmsearch) is applied to this profile.
To applied the --e-value-search to all profiles use the --no-cut-ga option."""},

               'i_evalue_sel': {'question': "Maximal independent e-value for Hmmer hits to be selected for systems detection.",
                                'validator': check_float,
                                'default': defaults.i_evalue_sel,
                                'explanation': ""},
               'coverage_profile': {'question': "Minimal profile coverage",
                                    'validator': check_float,
                                    'default': defaults.coverage_profile,
                                    'explanation':
"""Minimal profile coverage required for the hit alignment
with the profile to allow the hit selection for systems detection."""}
               }
    if not use_defaults:
        enter = ask("Do you want to enter Hmmer section?",
                    check_choice,
                    expected=["Y", "n"],
                    default="Y",
                    explanation="where to find hmmsearch, evalue, coverage, ...",
                    question_color=theme.EMPHASIZE + theme.SECTION
                    )
        use_defaults = enter == "n"
    set_section('hmmer', options, config, defaults, use_defaults=use_defaults)


def set_general_options(config, defaults, use_defaults=False):
    """
    Options for general section

    :param config: The config to setup
    :type config: :class:`ConfigParserWithComments` object
    :param defaults: the macsyfinder defaults values
    :type defaults: :class:`macsypy.config.MacsyDefaults` object
    :param bool use_defaults: If True do not ask any question use the defaults values
    """
    options = {'log_level': {'question': "The verbosity of the output",
                             'validator': check_choice,
                             'default': 'info',
                             'expected': ['debug', 'info', 'warning', 'error', 'critical'],
                             'explanation': ""},
               'worker': {'question': "Number of workers to be used by MacSyFinder.",
                          'validator': check_positive_int,
                          'default': defaults.worker,
                          'explanation':
"""In the case the user wants to run MacSyFinder in a multi-thread mode.
0 mean than one process by type of gene will be launch in parallel.."""},

               'mute': {'question': "Mute the log on stdout.",
                        'validator': check_bool,
                        'default': 'No',
                        'expected': ['Yes', 'No'],
                        'explanation': "Nothing is write in stdout, but MSF continue to log on macsyfinder.log"}
               }
    if not use_defaults:
        enter = ask("Do you want to enter general section?",
                    check_choice,
                    expected=["Y", "n"],
                    default="Y",
                    explanation="number of cpu used, verbosity, ...",
                    question_color=theme.EMPHASIZE + theme.SECTION
                    )
        use_defaults = enter == "n"
    set_section('general', options, config, defaults, use_defaults=use_defaults)


def set_score_options(config, defaults, use_defaults=False):
    """
    Options for scoring section

    :param config: The config to setup
    :type config: :class:`ConfigParserWithComments` object
    :param defaults: the macsyfinder defaults values
    :type defaults: :class:`macsypy.config.MacsyDefaults` object
    :param bool use_defaults: If True do not ask any question use the defaults values
    """
    options = {'mandatory_weight': {'question': "The weight of a mandatory component in cluster scoring.",
                                    'validator': check_float,
                                    'default': defaults.mandatory_weight,
                                    'explanation': ""},
               'accessory_weight': {'question': "The weight of a accessory component in cluster scoring.",
                                    'validator': check_float,
                                    'default': defaults.accessory_weight,
                                    'explanation': ""},
               'exchangeable_weight': {'question':
                                "The weight modifier for a component which code for exchangeable cluster scoring.",
                                       'validator': check_float,
                                       'default': defaults.exchangeable_weight,
                                       'explanation': ""},
               'redundancy_penalty': {'question':
                        "The weight modifier for cluster which bring a component already presents in an other one.",
                                      'validator': check_float,
                                      'default': defaults.redundancy_penalty,
                                      'explanation': ""},
               'out_of_cluster_weight': {'question': "The weight modifier for a hit which is not in a cluster",
                                         'validator': check_float,
                                         'default': defaults.out_of_cluster_weight,
                                         'explanation': """The hit is a
    - true loner (not in any cluster)
    - or multi-system (in a cluster but from an other system)"""}
               }
    if not use_defaults:
        enter = ask("Do you want to enter score section?",
                    check_choice,
                    expected=["Y", "n"],
                    default="Y",
                    explanation="The weights for mandatory, accessory, ...",
                    question_color=theme.EMPHASIZE + theme.SECTION
                    )
        use_defaults = enter == "n"
    set_section('score_opt', options, config, defaults, use_defaults=use_defaults)


def set_base_options(config, defaults, use_defaults=False):
    """
    Options for base section

    :param config: The config to setup
    :type config: :class:`ConfigParserWithComments` object
    :param defaults: the macsyfinder defaults values
    :type defaults: :class:`macsypy.config.MacsyDefaults` object
    :param bool use_defaults: If True do not ask any question use the defaults values
    """
    options = {'db_type': {'question': "The type sequence to analyze",
                           'validator': check_choice,
                           'default': str(defaults.db_type),
                           'expected': ['ordered_replicon', 'gembase', 'unordered', 'None'],
                           'explanation': ""},
               'replicon_topology': {'question': "The topology of replicon in dataset",
                                     'validator': check_choice,
                                     'default': defaults.replicon_topology,
                                     'expected': ['circular', 'linear'],
                                     'explanation': ""},

               'sequence_db': {'question': "The path to the sequence file.",
                               'validator': check_file,
                               'default': str(defaults.sequence_db),
                               'explanation': """By default macsyfinder will analyze this file.
But you can still specify another sequence file with --sequence-db option."""}
               }
    if not use_defaults:
        enter = ask("Do you want to enter in base section?",
                    check_choice,
                    expected=["Y", "n"],
                    default="Y",
                    explanation="Type of sequence to analyze, replicon topology, ...",
                    question_color=theme.EMPHASIZE + theme.SECTION
                    )
        use_defaults = enter == "n"
    set_section('base', options, config, defaults, use_defaults=use_defaults)


def prolog():
    """return the text displayed to the user when the configuration file is generated"""
    rep = f"""{theme.EMPHASIZE}Welcome to the MacSyFinder {msf_vers} configuration utility.{theme.RESET}

Please enter values for the following settings (just press Enter to
accept a default value, if one is given in brackets).
"""
    return rep


def epilog(path):
    """return the text to the user before to start the configuration"""
    rep = f"""A configuration file '{theme.EMPHASIZE}{path}{theme.RESET}' has been generated..
Place it in canonical location
 {theme.QUESTION}*{theme.RESET} in /etc/macsyfinder for system wide configuration {theme.WARN}(must named macsyfinder.conf){theme.RESET}
 {theme.QUESTION}*{theme.RESET} in <VIRTUALENV>/etc if you use a virtualenv {theme.WARN}(must named macsyfinder.conf){theme.RESET}
 {theme.QUESTION}*{theme.RESET} in ~/.macsyfinder for user wide configuration {theme.WARN}(must named macsyfinder.conf){theme.RESET}
 {theme.QUESTION}*{theme.RESET} where you run the analysis for local configuration {theme.WARN}(must named macsyfinder.conf){theme.RESET}
 {theme.QUESTION}*{theme.RESET} you can also put anywhere on the filesystems and use {theme.WARN}MACSY_CONF{theme.RESET} environment variable
   to indicate where to find it or specify it on the macsyfinder command line with option {theme.WARN}--cfg-file{theme.RESET}
   can be named as you want.

"""
    return rep


def serialize(config, path):
    """
    save the configuration on file

    :param config: the config to save
    :type config: :class:`ConfigParserWithComments` object
    :param str path: where to store the configuration
    """
    with open(path, 'w') as file:
        config.write(file)


def parse_args(args):
    """
    parse command line

    :param args: the command line arguments
    :type args: list of string
    :return:
    :rtype: :class:`argparse.Namespace` object
    """
    parser = argparse.ArgumentParser()
    theme_option = parser.add_mutually_exclusive_group()
    theme_option.add_argument("--no-color",
                              action='store_true',
                              default=False)
    theme_option.add_argument("--white-bg",
                              action='store_true',
                              default=False)
    theme_option.add_argument("--dark-bg",
                              action='store_true',
                              default=True)
    parser.add_argument("--defaults",
                        action='store_true',
                        default=False,
                        help="Do not ask questions. Create config file with default values.")
    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args=None) -> None:
    """
    The main entrypoint of the script

    :param args:
    """
    args = sys.argv[1:] if args is None else args
    parsed_args = parse_args(args)

    col_init()
    global theme
    if parsed_args.no_color:
        theme = Theme(ERROR=Style.BRIGHT,
                      WARN='',
                      SECTION='',
                      RESET=Style.RESET_ALL,
                      RETRY='',
                      QUESTION='',
                      EMPHASIZE='',
                      EXPLANATION='',
                      DEFAULT=''
                      )

    elif parsed_args.white_bg:
        theme = Theme(ERROR=Style.BRIGHT + Fore.RED,
                      WARN=Fore.LIGHTRED_EX,
                      RESET=Style.RESET_ALL,
                      RETRY=Style.BRIGHT + Fore.MAGENTA,
                      QUESTION=Style.BRIGHT,
                      EMPHASIZE=Style.BRIGHT,
                      EXPLANATION=Style.RESET_ALL,
                      DEFAULT=Style.BRIGHT + Fore.LIGHTBLACK_EX
                      )

    else:
        # parsed_args.dark_bg is always True
        # add in options only for coherence
        # and is the default
        pass

    config = ConfigParserWithComments()
    defaults = MacsyDefaults()
    conf_path = 'macsyfinder.conf'

    if os.path.exists(conf_path):
        go_on = ask(f"The '{conf_path}' already exists Overwrite/abort", check_choice,
                    expected=["O", "a"],
                    default="O",
                    question_color=theme.QUESTION)
        if go_on == "a":
            sys.exit(1)
    print(prolog())
    set_path_options(config, defaults, use_defaults=parsed_args.defaults)
    set_hmmer_options(config, defaults, use_defaults=parsed_args.defaults)
    set_score_options(config, defaults, use_defaults=parsed_args.defaults)
    set_general_options(config, defaults, use_defaults=parsed_args.defaults)
    set_base_options(config, defaults, use_defaults=parsed_args.defaults)
    serialize(config, conf_path)
    print(epilog(conf_path))


if __name__ == "__main__":
    main()
