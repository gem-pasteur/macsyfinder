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
import itertools

from operator import attrgetter  # To be used with "sorted"
from textwrap import dedent

import colorlog
_log = colorlog.getLogger('macsypy')

import macsypy
from macsypy.config import MacsyDefaults, Config
from macsypy.registries import ModelRegistry
from macsypy.definition_parser import DefinitionParser
from macsypy.search_genes import search_genes
from macsypy.database import Indexes, RepliconDB
#from macsypy.search_systems import search_systems
from macsypy.error import OptionError
from macsypy import cluster
from macsypy.hit import HitRegistry
from macsypy.system import match, System
from macsypy.utils import get_models_name_to_detect


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
                                     The first value must match to a model, the second to a number of components.\
                                     This option can be repeated several times:\
                                     \n \"--inter-gene-max-space T2SS 12 --inter-gene-max-space Flagellum 20\""
                                )
    system_options.add_argument("--min-mandatory-genes-required",
                                action='append',
                                nargs=2,
                                default=None,
                                help="The minimal number of mandatory genes required for model assessment.\
                                     The first value must correspond to a model name, the second value to an integer.\
                                     This option can be repeated several times:\n\
                                     \"--min-mandatory-genes-required T2SS 15 --min-mandatory-genes-required Flagellum 10\""
                                )
    system_options.add_argument("--min-genes-required",
                                action='append',
                                nargs=2,
                                default=None,
                                help="The minimal number of genes required for model assessment\
                                     (includes both 'mandatory' and 'accessory' components).\
                                     The first value must correspond to a model name, the second value to an integer.\
                                     This option can be repeated several times:\
                                     \n \"--min-genes-required T2SS 15 --min-genes-required Flagellum 10\""
                                )
    system_options.add_argument("--max-nb-genes",
                                action='append',
                                nargs=2,
                                default=None,
                                help="The maximal number of genes required for model assessment.\
                                     The first value must correspond to a model name, the second value to an integer.\
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
                               help='Maximal independent e-value for Hmmer hits to be selected for model detection.\
                                     (default = 0.001)')
    hmmer_options.add_argument('--coverage-profile',
                               action='store',
                               type=float,
                               default=None,
                               help='Minimal profile coverage required in the hit alignment to allow\
                                    the hit selection for model detection. (default = 0.5)')

    dir_options = parser.add_argument_group(title="Path options", description=None)
    dir_options.add_argument('--models-dir',
                             action='store',
                             default=None,
                             help='specify the path to the models if the models are not installed in the canonical place.\
                              It gather definitions (xml files) and hmm profiles in a specific\
                              structure. A directory with the name of the model with at least two directories\
                              "profiles" which contains all hmm profile for gene describe in definitions and\
                              "models" which contains either xml file of definitions or subdirectories\
                              to organize the model in subsystems.')
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

    parsed_args = parser.parse_args(args)
    return parsed_args


def main_search_systems(config, model_bank, gene_bank, profile_factory, logger):
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
    parser = DefinitionParser(config, model_bank, gene_bank, profile_factory)
    try:
        models_name_to_detect = get_models_name_to_detect(config.models(), registry)
    except KeyError as err:
        sys.exit("macsyfinder: {}".format(str(err).strip('"')))

    parser.parse(models_name_to_detect)

    logger.info("MacSyFinder's results will be stored in {0}".format(working_dir))
    logger.info("Analysis launched on {0} for model(s):".format(config.sequence_db()))

    for s in models_name_to_detect:
        logger.info("\t- {}".format(s))

    models_to_detect = [model_bank[model_fqn] for model_fqn in models_name_to_detect]
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
        # It's important to keep this sorting to have in last all_hits version
        # the hits with the same replicon_name and position sorted by score
        # the best score in first
        hits_by_replicon = {}
        for hit in all_hits:
            if hit.replicon_name in hits_by_replicon:
                hits_by_replicon[hit.replicon_name].append(hit)
            else:
                hits_by_replicon[hit.replicon_name] = [hit]

        for rep_name in hits_by_replicon:
            hits_by_replicon[rep_name].sort(key=attrgetter('score'), reverse=True)
            hits_by_replicon[rep_name].sort(key=attrgetter('position'))

        models_to_detect = sorted(models_to_detect, key=attrgetter('name'))
        rep_db = RepliconDB(config)
        hit_registry = HitRegistry()
        systems = []
        rejected_clusters = []
        for rep_name in hits_by_replicon:
            logger.info("Hits analysis for replicon {}".format(rep_name))
            rep_info = rep_db[rep_name]
            for model in models_to_detect:
                logger.info("Check model {}".format(model.fqn))
                hits_related_one_model = model.filter(hits_by_replicon[rep_name])
                logger.debug("{:#^80}".format(" hits related to {} ".format(model.name)))
                logger.debug("".join([str(h) for h in hits_related_one_model]))
                logger.debug("#" * 80)
                logger.info("Building clusters")
                clusters = cluster.build_clusters(hits_related_one_model, rep_info, model)
                logger.debug("{:#^80}".format("CLUSTERS"))
                logger.debug("\n".join([str(c) for c in clusters]))
                logger.debug("#" * 80)
                if model.multi_loci:
                    clusters_combination = [itertools.combinations(clusters, i) for i in range(1, len(clusters) + 1)]
                else:
                    clusters_combination = [clusters]
                logger.info("Searching systems")
                for one_combination_set in clusters_combination:
                    for one_clust_combination in one_combination_set:
                        res, _ = match(one_clust_combination, model, hit_registry)
                        if isinstance(res, System):
                            systems.append(res)
                        else:
                            rejected_clusters.append(res)

        system_filename = os.path.join(config.working_dir(), "macsyfinder.systems")
        systems_to_file(systems, system_filename)

        cluster_filename = os.path.join(config.working_dir(), "macsyfinder.rejected_cluster")
        rejected_clst_to_file(rejected_clusters, cluster_filename)
    else:
        logger.info("No hits found in this dataset.")


def systems_to_file(systems, system_filename):
    """
    print systems occurrences in a file

    :param systems: list of systems found
    :type systems: list of :class:`macsypy.system.System` objects
    :param str system_filename: The path of systems result file
    :return: None
    """
    with open(system_filename, "w") as sys_file:
        print("# macsyfinder {}".format(macsypy.__version__), file=sys_file)
        print("# {}".format(' '.join(sys.argv)), file=sys_file)
        print("# Systems found:\n", file=sys_file)
        for system in systems:
            print(system, file=sys_file)
            print("=" * 60, file=sys_file)


def rejected_clst_to_file(rejected_clusters, cluster_filename):
    """
    print rejected clusters in a file
    
    :param rejected_clusters: list of clusters which does not contitute a system
    :type rejected_clusters: list of :class:`macsypy.cluster.RejectedClusters` objects
    :param str cluster_filename: The path of rejected clusters result file
    :return: None
    """
    with open(cluster_filename, "w") as clst_file:
        print("# macsyfinder {}".format(macsypy.__version__), file=clst_file)
        print("# {}".format(' '.join(sys.argv)), file=clst_file)
        print("# Rejected clusters:\n", file=clst_file)

        for rej_clst in rejected_clusters:
            print(rej_clst, file=clst_file)
            print("=" * 60, file=clst_file)


def main(args=None, loglevel=None, models=None, genes=None, profiles=None):
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
        _log.info("command used: {}".format(' '.join(sys.argv)))
        if models is None:
            models = macsypy.model.ModelBank()
        if genes is None:
            genes = macsypy.gene.GeneBank()
        if profiles is None:
            profiles = macsypy.gene.ProfileFactory()
        main_search_systems(config, models, genes, profiles, logger)
    logger.debug("END")


if __name__ == "__main__":
    main()












