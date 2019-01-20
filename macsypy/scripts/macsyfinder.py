#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur, Paris.                                   #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


import sys
import os
import argparse
import logging
from operator import attrgetter  # To be used with "sorted"
from textwrap import dedent
import colorlog
_log = colorlog.getLogger('macsypy')

import macsypy
from macsypy.config import MacsyDefaults, Config
from macsypy.registries import ModelRegistry
from macsypy.system_parser import SystemParser
from macsypy.search_genes import search_genes
from macsypy.database import Indexes
from macsypy.search_systems import search_systems
from macsypy.system import system_bank
from macsypy.gene import gene_bank
from macsypy.error import OptionError



def get_models_name_to_detect(models, model_registry):
    """
    :param models: the list of models to detect as returned by config.models.
    :type models: list of tuple with the following structure:
                  [('model_1', ('def1, def2, ...)), ('model_2', ('def1', ...)), ...]
    :param model_registry: the models registry for this run.
    :type model_registry: :class:`macsypy.registries.ModelRegistry` object.
    :return: the models fully qualified name to launch a detection on.
    :rtype: list of string ['model_1/def1', 'model_1/def2', ..., 'model_2/def1', ...]
    :raise ValueError: if a model name provided in models is not in model_registry.
    """
    models_name_to_detect = []
    for group_of_defs in models:
        root = group_of_defs[0]
        definitions = group_of_defs[1]
        model_loc = model_registry[root.split('/')[0]]
        if 'all' in [d.lower() for d in definitions]:
            if root == model_loc.name:
                root = None
            def_loc = model_loc.get_all_definitions(root_def_name=root)
            models_name_to_detect.extend([d.fqn for d in def_loc])
        else:
            models_name_to_detect.extend([model_loc.get_definition('{}/{}'.format(root, one_def)).fqn
                                          for one_def in definitions])
    return models_name_to_detect


def get_version_message():
    """
    :return: the long description of the macsyfinder version
    :rtype: str
    """
    version = macsypy.__version__
    vers_msg = """Macsyfinder {0}
Python {1}

MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
See the COPYING file for details.

If you use this software please cite:
Abby SS, Néron B, Ménager H, Touchon M, Rocha EPC (2014)
MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems.
PLoS ONE 9(10): e110726. doi:10.1371/journal.pone.0110726""".format(version, sys.version)
    return vers_msg


def list_models(args):
    """
    :param args: The command line argument once parsed
    :type args: :class:`argparse.Namespace` object
    :return: a string representation of all models and submodels installed.
    :rtype: str
    """
    config = Config(MacsyDefaults(), args)
    registry = ModelRegistry(config)
    return str(registry)


def parse_args(args):
    """

    :param args:
    :return:
    """
    parser = argparse.ArgumentParser(
        epilog="For more details, visit the MacSyFinder website and see the MacSyFinder documentation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=dedent('''



         *            *               *                   *
    *           *               *   *   *  *    **                *   *
      **     *    *   *  *     *                    *               *
        __  __  *              ____ *        *  *  *    **     *
    || |  \/  | __ _  ___  || / ___| _   _  ||   ___ _         _        *
    || | |\/| |/ _` |/ __| || \___ \| | | | ||  | __(_)_ _  __| |___ _ _
    || | |  | | (_| | (__  ||  ___) | |_| | ||  | _|| | ' \/ _` / -_) '_|
    || |_|  |_|\__,_|\___| || |____/ \__, | ||  |_| |_|_||_\__,_\___|_|
               *             *       |___/         *                   *
     *      *   * *     *   **         *   *  *           *
      *      *         *        *    *              *
                 *                           *  *           *     *


    MacSyFinder - Detection of macromolecular systems in protein datasets 
    using systems modelling and similarity search.  
    '''))

    # , formatter_class=argparse.RawDescriptionHelpFormatter)
    # , formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", "--models",
                        action='append',
                        nargs='*',
                        default=None,
                        help='TODO bla-bla')

    genome_options = parser.add_argument_group(title="Input dataset options")
    genome_options.add_argument("--sequence-db",
                                action='store',
                                default=None,
                                help="Path to the sequence dataset in fasta format.")

    genome_options.add_argument("--db-type",
                                choices=['unordered_replicon', 'ordered_replicon', 'gembase', 'unordered'],
                                default=None,
                                help="The type of dataset to deal with. \"unordered_replicon\" corresponds\
                                     to a non-assembled genome,\"unordered\" to a metagenomic dataset,\
                                      \"ordered_replicon\" to an assembled genome, \
                                     and \"gembase\" to a set of replicons where sequence identifiers\
                                      follow this convention: \">RepliconName SequenceID\".")

    genome_options.add_argument("--replicon-topology",
                                choices=['linear', 'circular'],
                                default=None,
                                help="The topology of the replicons \
                                (this option is meaningful only if the db_type is 'ordered_replicon' or 'gembase'. ")

    genome_options.add_argument("--topology-file",
                                default=None,
                                help="Topology file path. The topology file allows to specify a topology \
                                     (linear or circular) for each replicon (this option is meaningful only if\
                                     the db_type is 'ordered_replicon' or 'gembase'.\
                                     A topology file is a tabular file with two columns: the 1st is the replicon name,\
                                     and the 2nd the corresponding topology:\n\"RepliconA\tlinear\" ")

    genome_options.add_argument("--idx",
                                action='store_true',
                                default=False,
                                help="Forces to build the indexes for the sequence dataset even \
                                     if they were previously computed and present at the dataset location (default = False)"
                                )

    system_options = parser.add_argument_group(title="Systems detection options")
    system_options.add_argument("--inter-gene-max-space",
                                action='append',
                                nargs=2,
                                default=None,
                                help="Co-localization criterion: maximum number of components non-matched by a\
                                     profile allowed between two matched components for them to be considered contiguous.\
                                     Option only meaningful for 'ordered' datasets.\
                                     The first value must match to a system, the second to a number of components.\
                                     This option can be repeated several times:\
                                     \n \"--inter-gene-max-space T2SS 12 --inter-gene-max-space Flagellum 20\""
                                )
    system_options.add_argument("--min-mandatory-genes-required",
                                action='append',
                                nargs=2,
                                default=None,
                                help="The minimal number of mandatory genes required for system assessment.\
                                     The first value must correspond to a system name, the second value to an integer.\
                                     This option can be repeated several times:\n\
                                     \"--min-mandatory-genes-required T2SS 15 --min-mandatory-genes-required Flagellum 10\""
                                )
    system_options.add_argument("--min-genes-required",
                                action='append',
                                nargs=2,
                                default=None,
                                help="The minimal number of genes required for system assessment\
                                     (includes both 'mandatory' and 'accessory' components).\
                                     The first value must correspond to a system name, the second value to an integer.\
                                     This option can be repeated several times:\
                                     \n \"--min-genes-required T2SS 15 --min-genes-required Flagellum 10\""
                                )
    system_options.add_argument("--max-nb-genes",
                                action='append',
                                nargs=2,
                                default=None,
                                help="The maximal number of genes required for system assessment.\
                                     The first value must correspond to a system name, the second value to an integer.\
                                     This option can be repeated several times:\
                                     \n \"--max-nb-genes T2SS 5 --max-nb-genes Flagellum 10"
                                )
    system_options.add_argument("--multi-loci",
                                action='store',
                                default=None,
                                help="Allow the storage of multi-loci systems for the specified systems.\
                                The systems are specified as a comma separated list \
                                (--multi-loci sys1,sys2) default is False"
                                )

    hmmer_options = parser.add_argument_group(title="Options for Hmmer execution and hits filtering")
    hmmer_options.add_argument('--hmmer',
                               action='store',
                               default=None,
                               help='Path to the Hmmer program.')
    hmmer_options.add_argument('--e-value-search',
                               action='store',
                               type=float,
                               default=None,
                               help='Maximal e-value for hits to be reported during Hmmer search. (default = 1)')
    hmmer_options.add_argument('--i-evalue-sel',
                               action='store',
                               type=float,
                               default=None,
                               help='Maximal independent e-value for Hmmer hits to be selected for system detection.\
                                     (default = 0.001)')
    hmmer_options.add_argument('--coverage-profile',
                               action='store',
                               type=float,
                               default=None,
                               help='Minimal profile coverage required in the hit alignment to allow\
                                    the hit selection for system detection. (default = 0.5)')

    dir_options = parser.add_argument_group(title="Path options", description=None)
    dir_options.add_argument('--models-dir',
                             action='store',
                             default=None,
                             help='specify the path to the models if the models are not installed in the canonical place.\
                              It gather definitions (xml files) and hmm profiles in a specific\
                              structure. A directory with the name of the system with at least two directories\
                              "profiles" which contains all hmm profile for gene describe in definitions and\
                              "models" which contains either xml file of definitions or subdirectories\
                              to organize the system in subsystems.')
    dir_options.add_argument('-o', '--out-dir',
                             action='store',
                             default=None,
                             help='Path to the directory where to store results.\
                             if out-dir is specified res-search-dir will be ignored.')
    dir_options.add_argument('--res-search-suffix',
                             action='store',
                             default=None,
                             help='The suffix to give to Hmmer raw output files.')
    dir_options.add_argument('--res-extract-suffix',
                             action='store',
                             default=None,
                             help='The suffix to give to filtered hits output files.')
    dir_options.add_argument('--profile-suffix',
                             action='store',
                             default=None,
                             help="The suffix of profile files. For each 'Gene' element, the corresponding profile is \
                                  searched in the 'profile_dir', in a file which name is based on the \
                                  Gene name + the profile suffix.\
                                  For instance, if the Gene is named 'gspG' and the suffix is '.hmm3',\
                                  then the profile should be placed at the specified location and be named 'gspG.hmm3'")

    general_options = parser.add_argument_group(title="General options", description=None)
    general_options.add_argument("-w", "--worker",
                                 action='store',
                                 type=int,
                                 default=None,
                                 help="Number of workers to be used by MacSyFinder.\
                                      In the case the user wants to run MacSyFinder in a multi-thread mode.\
                                      (0 mean all cores will be used, default 1)")
    general_options.add_argument("-v", "--verbosity",
                                 action="count",
                                 default=0,
                                 help="Increases the verbosity level. There are 4 levels:\
                                       Error messages (default), Warning (-v), Info (-vv) and Debug.(-vvv)")
    general_options.add_argument("--mute",
                                 action="store_true",
                                 default=False,
                                 help="mute the log on stdout."
                                      " (continue to log on macsyfinder.log)")
    general_options.add_argument("--version",
                                 action="version",
                                 version=get_version_message()),
    general_options.add_argument("-l", "--list-models",
                                 action="store_true",
                                 default=False,
                                 help="display the all models installed in generic location and quit.")
    general_options.add_argument("--cfg-file",
                                 action='store',
                                 help="Path to a MacSyFinder configuration file to be used.")
    general_options.add_argument("--previous-run",
                                 action='store',
                                 default=None,
                                 help="""Path to a previous MacSyFinder run directory.
                                         It allows to skip the Hmmer search step on same dataset,
                                         as it uses previous run results and thus parameters regarding Hmmer detection.
                                         The configuration file from this previous run will be used.
                                         (conflict with options  --config, --sequence-db, --profile-suffix,
                                         --res-extract-suffix, --e-value-res, --db-type, --hmmer)""")
    general_options.add_argument("--relative-path",
                                 action='store_true',
                                 default=False,
                                 help=argparse.SUPPRESS)
    # 'relative-path' option help message (currently hidden)
    # Use relative paths instead of absolute paths. This option is used
    # by developers to generate portable data set, as for example test
    # data set, which are used on many different machines (using previous-run option).

    parsed_args = parser.parse_args()
    return parsed_args


def main_search_systems(config, logger):
    """

    :param parsed_args: the command line arguments
    :type parsed_args: a :class:`argparse.Namespace` object
    :param logger:
    :return:
    """
    working_dir = config.working_dir()
    config.save(path_or_buf=os.path.join(working_dir, config.cfg_name))

    registry = ModelRegistry(config)
    # build indexes
    idx = Indexes(config)
    idx.build(force=config.idx)

    # create models
    parser = SystemParser(config, system_bank, gene_bank)
    try:
        models_name_to_detect = get_models_name_to_detect(config.models(), registry)
    except KeyError as err:
        sys.exit("macsyfinder: {}".format(str(err).strip('"')))

    parser.parse(models_name_to_detect)

    logger.info("MacSyFinder's results will be stored in {0}".format(working_dir))
    logger.info("Analysis launched on {0} for system(s):".format(config.sequence_db()))

    for s in models_name_to_detect:
        logger.info("\t- {}".format(s))

    models_to_detect = [system_bank[model_fqn] for model_fqn in models_name_to_detect]
    all_genes = []
    for system in models_to_detect:
        genes = system.mandatory_genes + system.accessory_genes + system.forbidden_genes
        # Exchangeable homologs/analogs are also added cause they can "replace" an important gene...
        ex_genes = []

        for g in genes:
            if g.exchangeable:
                h_s = g.get_homologs()
                ex_genes += h_s
                a_s = g.get_analogs()
                ex_genes += a_s
        all_genes += (genes + ex_genes)
    #############################################
    # this part of code is executed in parallel
    #############################################
    try:
        all_reports = search_genes(all_genes, config)
    except Exception as err:
        sys.exit(str(err))
    #############################################
    # end of parallel code
    #############################################
    all_hits = [hit for subl in [report.hits for report in all_reports] for hit in subl]

    if len(all_hits) > 0:
        # It's important to keep this sorting to have in last  all_hits version
        # the hits with the same replicon_name and position sorted by score
        # the best score in first
        all_hits = sorted(all_hits, key=attrgetter('score'), reverse=True)
        all_hits = sorted(all_hits, key=attrgetter('replicon_name', 'position'))
        models_to_detect = sorted(models_to_detect, key=attrgetter('name'))
        search_systems(all_hits, models_to_detect, config)
    else:
        logger.info("No hits found in this dataset.")


def main(args=None, loglevel=None):
    """
    main entry point to integron_finder

    :param str args: the arguments passed on the command line
    :param loglevel: the output verbosity
    :type loglevel: a positive int or a string among 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    args = sys.argv[1:] if args is None else args
    parsed_args = parse_args(args)

    defaults = MacsyDefaults()
    config = Config(defaults, parsed_args)

    ###########################
    # creation of working dir
    ###########################

    working_dir = config.working_dir()
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    else:
        if os.path.isdir(working_dir):
            if os.listdir(working_dir):
                raise ValueError("'{}' already exists and is not a empty".format(working_dir))
        else:
            raise ValueError("'{}' already exists and is not a directory".format(working_dir))

    ################
    # init loggers #
    ################
    macsypy.init_logger(log_file=os.path.join(config.working_dir(), config.log_file()),
                        out=not config.mute())
    if not loglevel:
        # logs are specify from args options
        macsypy.logger_set_level(level=config.log_level())
    else:
        # used by unit tests to mute or unmute logs
        macsypy.logger_set_level(level=loglevel)

    logger = logging.getLogger('macsypy.macsyfinder')


    if parsed_args.list_models:
        print(list_models(parsed_args), file=sys.stdout)
        sys.exit(0)
    else:
        if not parsed_args.previous_run and not parsed_args.models:
            raise OptionError("argument --models is required.")
        if not parsed_args.previous_run and not parsed_args.sequence_db:
            raise OptionError("argument --sequence-db is required.")
        if not parsed_args.previous_run and not parsed_args.db_type:
            raise OptionError("argument --db-type is required.")

        main_search_systems(config, logger)
    logger.debug("END")

if __name__ == "__main__":

    main()












