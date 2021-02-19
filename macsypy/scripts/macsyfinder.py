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
from macsypy.registries import ModelRegistry, scan_models_dir
from macsypy.definition_parser import DefinitionParser
from macsypy.search_genes import search_genes
from macsypy.database import Indexes, RepliconDB
from macsypy.error import MacsypyError, OptionError
from macsypy import cluster
from macsypy.hit import get_best_hits, HitWeight
from macsypy.system import OrderedMatchMaker, UnorderedMatchMaker, System, LikelySystem, HitSystemTracker
from macsypy.utils import get_def_to_detect
from macsypy.profile import ProfileFactory
from macsypy.model import ModelBank
from macsypy.gene import GeneBank
from macsypy.solution import find_best_solutions
from macsypy.serialization import TxtSystemSerializer, TxtLikelySystemSerializer, TxtUnikelySystemSerializer, \
    TsvSystemSerializer, TsvSolutionSerializer, TsvLikelySystemSerializer


def get_version_message():
    """
    :return: the long description of the macsyfinder version
    :rtype: str
    """
    version = macsypy.__version__
    py_vers = sys.version.replace('\n', ' ')
    vers_msg = f"""Macsyfinder {version}
using:
- Python {py_vers}
- NetworkX {macsypy.solution.nx.__version__}

MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
See the COPYING file for details.

If you use this software please cite:
{macsypy.__citation__}
and don't forget to cite models used:
macsydata cite <model>
"""
    return vers_msg


def list_models(args):
    """
    :param args: The command line argument once parsed
    :type args: :class:`argparse.Namespace` object
    :return: a string representation of all models and submodels installed.
    :rtype: str
    """
    defaults = MacsyDefaults()
    config = Config(defaults, args)
    system_model_dir = config.models_dir()
    if args.models_dir:
        model_dirs = (system_model_dir,)
    else:
        user_model_dir = os.path.join(os.path.expanduser('~'), '.macsyfinder', 'data')
        model_dirs = (system_model_dir, user_model_dir) if os.path.exists(user_model_dir) else (system_model_dir,)
    registry = ModelRegistry()
    for model_dir in model_dirs:
        try:
            for model_loc in scan_models_dir(model_dir, profile_suffix=config.profile_suffix):
                registry.add(model_loc)
        except PermissionError as err:
            _log.warning(f"{model_dir} is not readable: {err} : skip it.")
    return str(registry)


def parse_args(args):
    """

    :param args: The arguments provided on the command line
    :type args: List of strings [without the program name]
    :return: The arguments parsed
    :rtype: :class:`aprgparse.Namespace` object.
    """
    parser = argparse.ArgumentParser(
        epilog="For more details, visit the MacSyFinder website and see the MacSyFinder documentation.",
        # formatter_class=ArgumentDefaultsHelpRawTextFormatter,
        formatter_class=argparse.RawTextHelpFormatter,
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


    MacSyFinder (MSF) - Detection of macromolecular systems in protein datasets 
    using systems modelling and similarity search.  
    '''))

    msf_def = MacsyDefaults()
    # , formatter_class=argparse.RawDescriptionHelpFormatter)
    # , formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-m", "--models",
                        nargs='*',
                        default=None,
                        help="""The models to search.
The first element must be the name of family models, followed by the name of the models to search.
If the name 'all' is in the list of models, all models from the family will be searched.
'--models TXSS Flagellum T2SS' 
          means MSF will search for the models TXSS/Flagellum and TXSS/T2SS
'--models TXSS all' 
          means MSF will search for all models found in the model package TXSS
'--models CRISPRcas/subtyping all' 
          means MSF will search for all models described in the CRISPRCas/subtyping subfamily.
(required unless --previous-run is set)
""")

    genome_options = parser.add_argument_group(title="Input dataset options")
    genome_options.add_argument("--sequence-db",
                                action='store',
                                default=None,
                                help="""Path to the sequence dataset in fasta format.
(required unless --previous-run is set)
""")

    genome_options.add_argument("--db-type",
                                choices=['ordered_replicon', 'gembase', 'unordered'],
                                default=None,
                                help='''The type of dataset to deal with. 
"unordered" corresponds to a non-assembled genome or set of unassembled genes,
"ordered_replicon" to an assembled genome,
"gembase" to a set of replicons where sequence identifiers
	follow this convention: ">RepliconName_SequenceID".
(required unless --previous-run is set)
''')

    genome_options.add_argument("--replicon-topology",
                                choices=['linear', 'circular'],
                                default=None,
                                help=f"""The topology of the replicons
(this option is meaningful only if the db_type is
'ordered_replicon' or 'gembase'.)
(default: {msf_def['e_value_search']})
""")

    genome_options.add_argument("--topology-file",
                                default=None,
                                help="""Topology file path. The topology file allows to specify a topology
(linear or circular) for each replicon (this option is meaningful only if the db_type is 
'ordered_replicon' or 'gembase'.
A topology file is a tabular file with two columns:
	the 1st is the replicon name, and the 2nd the corresponding topology:
	\"RepliconA\tlinear\"
""")

    genome_options.add_argument("--idx",
                                action='store_true',
                                default=False,
                                help=f"""Forces to build the indexes for the sequence dataset even
if they were previously computed and present at the dataset location.
(default: {msf_def['idx']})"""
                                )

    system_options = parser.add_argument_group(title="Systems detection options")
    system_options.add_argument("--inter-gene-max-space",
                                action='append',
                                nargs=2,
                                default=None,
                                help="""Co-localization criterion: maximum number of components non-matched by a
	profile allowed between two matched components for them to be considered contiguous.
Option only meaningful for 'ordered' datasets.
The first value must name a model, the second a number of components.
This option can be repeated several times:
    "--inter-gene-max-space TXSS/T2SS 12 --inter-gene-max-space TXSS/Flagellum 20
"""
                                )
    system_options.add_argument("--min-mandatory-genes-required",
                                action='append',
                                nargs=2,
                                default=None,
                                help="""The minimal number of mandatory genes required for model assessment.
The first value must correspond to a model fully qualified name, the second value to an integer.
This option can be repeated several times:
    "--min-mandatory-genes-required TXSS/T2SS 15 --min-mandatory-genes-required TXSS/Flagellum 10"
"""
                                )
    system_options.add_argument("--min-genes-required",
                                action='append',
                                nargs=2,
                                default=None,
                                help="""The minimal number of genes required for model assessment 
(includes both 'mandatory' and 'accessory' components).
The first value must correspond to a model fully qualified name, the second value to an integer.
This option can be repeated several times:
    "--min-genes-required TXSS/T2SS 15 --min-genes-required TXSS/Flagellum 10
"""
                                )
    system_options.add_argument("--max-nb-genes",
                                action='append',
                                nargs=2,
                                default=None,
                                help="""The maximal number of genes to consider a system as full.
The first value must correspond to a model name, the second value to an integer.
This option can be repeated several times:
    "--max-nb-genes TXSS/T2SS 5 --max-nb-genes TXSS/Flagellum 10"
"""
                                )
    system_options.add_argument("--multi-loci",
                                action='store',
                                default=None,
                                help="""Specifies if the system can be detected as a 'scattered' (or multiple-loci-encoded) system.
The models are specified as a comma separated list of fully qualified name(s)
    "--multi-loci model_familyA/model_1,model_familyB/model_2"
""")
    hmmer_options = parser.add_argument_group(title="Options for Hmmer execution and hits filtering")
    hmmer_options.add_argument('--hmmer',
                               action='store',
                               default=None,
                               help=f"""Path to the hmmsearch program.
If not specified, rely on the environment variable PATH
(default: {msf_def['hmmer']})""")
    hmmer_options.add_argument('--e-value-search',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""Maximal e-value for hits to be reported during hmmsearch search.
By default MSF set per profile threshold for hmmsearch run (hmmsearch --cut_ga option) 
for profiles containing the GA bit score threshold.
If a profile does not contains the GA bit score the --e-value-search (-E in hmmsearch) is applied to this profile.
To applied the --e-value-search to all profiles use the --no-cut-ga option. 
(default: {msf_def['e_value_search']})
""")
    hmmer_options.add_argument('--no-cut-ga',
                               action='store_true',
                               default=False,
                               help=f"""By default the Mf try to applied a threshold per profile by using the 
hmmer -cut-ga option. This is possible only if the Ga bit score is present in the profile otherwise MF switch to use the
the --e-value-search (-E in hmmsearch). 
If this option is set the --e-value-search option is used for all profiles regardless 
the presence of the a GA bit score in the profiles.
(default: {msf_def['no_cut_ga']})""")

    hmmer_options.add_argument('--i-evalue-sel',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""Maximal independent e-value for Hmmer hits to be selected for systems detection.
(default:{msf_def['i_evalue_sel']})""")
    hmmer_options.add_argument('--coverage-profile',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""Minimal profile coverage required for the hit alignment  with the profile to allow
the hit selection for systems detection. 
(default: {msf_def['coverage_profile']})""")
    score_options = parser.add_argument_group(title="Score options",
                                              description="Options for cluster and systems scoring")
    score_options.add_argument('--mandatory-weight',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""the weight of a mandatory component in cluster scoring
(default:{msf_def['mandatory_weight']})""")
    score_options.add_argument('--accessory-weight',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""the weight of a mandatory component in cluster scoring
(default:{msf_def['accessory_weight']})""")

    # the weight of a mandatory component in cluster scoring
    # (default:{msf_def['neutral_weight']})
    score_options.add_argument('--neutral-weight',
                               action='store',
                               type=float,
                               default=None,
                               help=argparse.SUPPRESS)

    # the weight modifier for a component which code for itself cluster scoring
    # (default:{msf_def['itself_weight']})"""
    score_options.add_argument('--itself-weight',
                               action='store',
                               type=float,
                               default=None,
                               help=argparse.SUPPRESS)

    score_options.add_argument('--exchangeable-weight',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""the weight modifier for a component which code for exchangeable cluster scoring
    (default:{msf_def['exchangeable_weight']})""")
    score_options.add_argument('--redundancy-penalty',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""the weight modifier for cluster which bring a component already presents in other 
clusters (default:{msf_def['redundancy_penalty']})""")
    score_options.add_argument('--loner-multi-system-weight',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""the weight modifier for a hit which is a loner and multi-system 
(default:{msf_def['loner_multi_system_weight']})""")

    dir_options = parser.add_argument_group(title="Path options", description=None)
    dir_options.add_argument('--models-dir',
                             action='store',
                             default=None,
                             help="""Specifies the path to the models if the models are not installed in the canonical place.
It gathers definitions (xml files) and HMM profiles arranged in a specific
file structure. A directory with the name of the model with at least two directories
	'profiles' - which contains HMM profiles for each gene components described in the systems' models
	'models' - which contains either the XML files of models' definitions or subdirectories
to organize the models in subsystems.""")
    dir_options.add_argument('-o', '--out-dir',
                             action='store',
                             default=None,
                             help="""Path to the directory where to store output results.
if out-dir is specified, res-search-dir will be ignored.""")
    dir_options.add_argument('--res-search-suffix',
                             action='store',
                             default=None,
                             help="The suffix to give to Hmmer raw output files. "
                                  f"(default: {msf_def['res_search_suffix']})")
    dir_options.add_argument('--res-extract-suffix',
                             action='store',
                             default=None,
                             help="The suffix to give to filtered hits output files. "
                                  f"(default: {msf_def['res_extract_suffix']})")
    dir_options.add_argument('--profile-suffix',
                             action='store',
                             default=None,
                             help=f"""The suffix of profile files. For each 'Gene' element, the corresponding profile is
searched in the 'profile_dir', in a file which name is based on the
Gene name + the profile suffix.
For instance, if the Gene is named 'gspG' and the suffix is '.hmm3',
then the profile should be placed at the specified location 
under the name 'gspG.hmm3'
(default: {msf_def['profile_suffix']})"""
                             )

    general_options = parser.add_argument_group(title="General options", description=None)
    general_options.add_argument("-w", "--worker",
                                 action='store',
                                 type=int,
                                 default=None,
                                 help=f"""Number of workers to be used by MacSyFinder.
In the case the user wants to run MacSyFinder in a multi-thread mode.
(0 mean all cores will be used).
(default: {msf_def['worker']})"""
                                 )
    general_options.add_argument("-v", "--verbosity",
                                 action="count",
                                 default=0,
                                 help="""Increases the verbosity level. There are 4 levels:
Error messages (default), Warning (-v), Info (-vv) and Debug.(-vvv)""")
    general_options.add_argument("--mute",
                                 action="store_true",
                                 default=False,
                                 help=f"""Mute the log on stdout.
(continue to log on macsyfinder.log)
(default: {msf_def['mute']})""")
    general_options.add_argument("--version",
                                 action="version",
                                 version=get_version_message())
    general_options.add_argument("-l", "--list-models",
                                 action="store_true",
                                 default=False,
                                 help="Displays all models installed at generic location and quit.")
    general_options.add_argument("--cfg-file",
                                 action='store',
                                 help="Path to a MacSyFinder configuration file to be used. (conflict with --previous-run)")
    general_options.add_argument("--previous-run",
                                 action='store',
                                 default=None,
                                 help="""Path to a previous MacSyFinder run directory.
It allows to skip the Hmmer search step on a same dataset,
as it uses previous run results and thus parameters regarding Hmmer detection.
The configuration file from this previous run will be used.
Conflicts with options:  
    --cfg-file, --sequence-db, --profile-suffix, --res-extract-suffix, --e-value-res, --db-type, --hmmer""")
    general_options.add_argument("--relative-path",
                                 action='store_true',
                                 default=False,
                                 help=argparse.SUPPRESS)
    # 'relative-path' option help message (currently hidden)
    # Use relative paths instead of absolute paths. This option is used
    # by developers to generate portable data set, as for example test
    # data set, which are used on many different machines (using previous-run option).

    parsed_args = parser.parse_args(args)
    if parsed_args.cfg_file and parsed_args.previous_run:
        # argparse does not allow to have mutually exclusive option  in a argument group
        # I prefer to have these 2 options in general options group
        # so I mimic the exclusive_group behavior
        parser.print_usage()
        print("macsyfinder: error: argument --previous-run: not allowed with argument --cfg-file")
        sys.exit(2)
    return parser, parsed_args


def search_systems(config, model_bank, gene_bank, profile_factory, logger):
    """
    Do the job, this function is the orchestrator of all the macsyfinder mechanics
    at the end several files are produced containing the results

      - macsyfinder.conf: The set of variables used to runt this job
      - macsyfinder.systems: The list of the potential systems
      - macsyfinder.rejected_cluster: The list of all clusters and clustrs combination
                                      which has been rejected and the reason
      - macsyfinder.log: the copy of the standard output

    :param config: The MacSyFinder Configuration
    :type config: :class:`macsypy.config.Config` object
    :param model_bank: The bank populated with the available models
    :type model_bank: :class:`macsypy.model.ModelBank` object
    :param gene_bank: the bank containing all genes
    :type gene_bank: :class:`macsypy.gene.GeneBank` object
    :param profile_factory: The profile factory
    :type profile_factory: :class:`macsypy.gene.ProfileFactory`
    :param logger: The logger use to display information to the user.
                   It must be initialized. see :func:`macsypy.init_logger`
    :type logger: :class:`colorlog.Logger` object
    :return: the systems and rejected clusters found
    :rtype: ([:class:`macsypy.system.System`, ...], [:class:`macsypy.cluster.RejectedCluster`, ...])
    """
    working_dir = config.working_dir()
    config.save(path_or_buf=os.path.join(working_dir, config.cfg_name))
    registry = ModelRegistry()
    models_loc_available = scan_models_dir(config.models_dir(),
                                           profile_suffix=config.profile_suffix(),
                                           relative_path=config.relative_path())
    for model_loc in models_loc_available:
        registry.add(model_loc)
    # build indexes
    idx = Indexes(config)
    idx.build(force=config.idx)

    # create models
    parser = DefinitionParser(config, model_bank, gene_bank, registry, profile_factory)
    try:
        models_def_to_detect = get_def_to_detect(config.models(), registry)
    except KeyError as err:
        sys.exit(f"macsyfinder: {err}")

    parser.parse(models_def_to_detect)

    logger.info(f"MacSyFinder's results will be stored in working_dir{working_dir}")
    logger.info(f"Analysis launched on {config.sequence_db()} for model(s):")

    for m in models_def_to_detect:
        logger.info(f"\t- {m.fqn}")

    models_to_detect = [model_bank[model_loc.fqn] for model_loc in models_def_to_detect]
    all_genes = []
    for model in models_to_detect:
        genes = model.mandatory_genes + model.accessory_genes + model.neutral_genes + model.forbidden_genes
        # Exchangeable (formerly homologs/analogs) are also added because they can "replace" an important gene...
        ex_genes = []

        for g in genes:
            ex_genes += g.exchangeables
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
            hits_by_replicon[rep_name] = get_best_hits(hits_by_replicon[rep_name], key='score')
            hits_by_replicon[rep_name].sort(key=attrgetter('position'))

        models_to_detect = sorted(models_to_detect, key=attrgetter('name'))
        db_type = config.db_type()
        if db_type in ('ordered_replicon', 'gembase'):
            systems, rejected_clusters = _search_in_ordered_replicon(hits_by_replicon, models_to_detect,
                                                                     config, logger)
            return systems, rejected_clusters
        elif db_type == "unordered":
            likely_systems, rejected_hits = _search_in_unordered_replicon(hits_by_replicon, models_to_detect,
                                                                          logger)
            return likely_systems, rejected_hits
        else:
            assert False, f"dbtype have an invalid value {db_type}"
    else:
        # No hits detected
        return [], []


def _search_in_ordered_replicon(hits_by_replicon, models_to_detect, config, logger):
    systems = []
    rejected_clusters = []
    rep_db = RepliconDB(config)
    for rep_name in hits_by_replicon:
        logger.info("\n{:#^60}".format(f" Hits analysis for replicon {rep_name} "))
        rep_info = rep_db[rep_name]
        for model in models_to_detect:
            logger.info(f"Check model {model.fqn}")
            hits_related_one_model = model.filter(hits_by_replicon[rep_name])
            logger.debug("{:#^80}".format(" hits related to {} ".format(model.name)))
            logger.debug("".join([str(h) for h in hits_related_one_model]))
            logger.debug("#" * 80)
            logger.info("Building clusters")
            hit_weights = HitWeight(**config.hit_weights())
            clusters = cluster.build_clusters(hits_related_one_model, rep_info, model, hit_weights)
            logger.debug("{:#^80}".format("CLUSTERS"))
            logger.debug("\n" + "\n".join([str(c) for c in clusters]))
            logger.debug("#" * 80)
            logger.info("Searching systems")
            if model.multi_loci:
                # The loners are already in clusters lists with their context
                # so they are take in account
                clusters_combination = [itertools.combinations(clusters, i) for i in range(1, len(clusters) + 1)]
            else:
                # we must add loners manually
                # but only if the cluster does not already contains them
                loners = cluster.get_loners(hits_related_one_model, model, hit_weights)
                clusters_combination = []
                for one_cluster in clusters:
                    one_clust_combination = [one_cluster]
                    filtered_loners = cluster.filter_loners(one_cluster, loners)
                    one_clust_combination.extend(filtered_loners)
                    clusters_combination.append([one_clust_combination])

            for one_combination_set in clusters_combination:
                for one_clust_combination in one_combination_set:
                    ordered_matcher = OrderedMatchMaker(model, redundancy_penalty=config.redundancy_penalty())
                    res = ordered_matcher.match(one_clust_combination)
                    if isinstance(res, System):
                        systems.append(res)
                    else:
                        rejected_clusters.append(res)
    if systems:
        systems.sort(key=lambda syst: (syst.replicon_name, syst.position[0], syst.model.fqn, - syst.score))
    return systems, rejected_clusters


def _search_in_unordered_replicon(hits_by_replicon, models_to_detect, logger):
    likely_systems = []
    rejected_hits = []
    for rep_name in hits_by_replicon:
        logger.info("\n{:#^60}".format(f" Hits analysis for replicon {rep_name} "))
        for model in models_to_detect:
            logger.info(f"Check model {model.fqn}")
            hits_related_one_model = model.filter(hits_by_replicon[rep_name])
            logger.debug("{:#^80}".format(" hits related to {} ".format(model.name)))
            logger.debug("".join([str(h) for h in hits_related_one_model]))
            logger.debug("#" * 80)
            logger.info("Searching systems")
            hits_related_one_model = model.filter(hits_by_replicon[rep_name])
            if hits_related_one_model:
                unordered_matcher = UnorderedMatchMaker(model)
                res = unordered_matcher.match(hits_related_one_model)
                if isinstance(res, LikelySystem):
                    likely_systems.append(res)
                else:
                    rejected_hits.append(res)
            else:
                logger.info(f"No hits found for model {model.fqn}")
    if likely_systems:
        likely_systems.sort(key=lambda syst: (syst.replicon_name, syst.position[0], syst.model.fqn))
    return likely_systems, rejected_hits


def _outfile_header():
    header = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}"""
    return header


def systems_to_tsv(systems, hit_system_tracker, sys_file):
    """
    print systems occurrences in a file in tabulated  format

    :param systems: list of systems found
    :type systems: list of :class:`macsypy.system.System` objects
    :param hit_system_tracker: a filled HitSystemTracker.
    :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :return: None
    """
    print(_outfile_header(), file=sys_file)
    if systems:
        print("# Systems found:", file=sys_file)
        print(TsvSystemSerializer.header, file=sys_file)
        for system in systems:
            sys_serializer = TsvSystemSerializer()
            print(sys_serializer.serialize(system, hit_system_tracker), file=sys_file)
    else:
        print("# No Systems found", file=sys_file)


def systems_to_txt(systems, hit_system_tracker, sys_file):
    """
    print systems occurrences in a file in human readable format

    :param systems: list of systems found
    :type systems: list of :class:`macsypy.system.System` objects
    :param hit_system_tracker: a filled HitSystemTracker.
    :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :return: None
    """

    print(_outfile_header(), file=sys_file)
    if systems:
        print("# Systems found:\n", file=sys_file)
        for system in systems:
            sys_serializer = TxtSystemSerializer()
            print(sys_serializer.serialize(system, hit_system_tracker), file=sys_file)
            print("=" * 60, file=sys_file)
    else:
        print("# No Systems found", file=sys_file)


def solutions_to_tsv(solutions, hit_system_tracker, sys_file):
    """
    print solution in a file in tabulated format
    A solution is a set of systems which represents an optimal combination of
    systems to maximize the score.

    :param solutions: list of systems found
    :type solutions: list of list of :class:`macsypy.system.System` objects
    :param hit_system_tracker: a filled HitSystemTracker.
    :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :return: None
    """
    print(_outfile_header(), file=sys_file)
    if solutions:
        sol_serializer = TsvSolutionSerializer()
        print("# Systems found:", file=sys_file)
        print(sol_serializer.header, file=sys_file)
        for sol_id, solution in enumerate(solutions, 1):
            solution.sort(key=lambda syst: (syst.replicon_name, syst.position[0], syst.model.fqn, - syst.score))
            print(sol_serializer.serialize(solution, sol_id, hit_system_tracker), file=sys_file, end='')
    else:
        print("# No Systems found", file=sys_file)


def rejected_clst_to_txt(rejected_clusters, clst_file):
    """
    print rejected clusters in a file

    :param rejected_clusters: list of clusters which does not contitute a system
    :type rejected_clusters: list of :class:`macsypy.cluster.RejectedClusters` objects
    :param clst_file: The file where to write down the rejected clusters
    :type clst_file: file object
    :return: None
    """
    print(_outfile_header(), file=clst_file)
    if rejected_clusters:
        print("# Rejected clusters:\n", file=clst_file)
        for rej_clst in rejected_clusters:
            print(rej_clst, file=clst_file, end='')
            print("=" * 60, file=clst_file)
    else:
        print("# No Rejected clusters", file=clst_file)


def likely_systems_to_txt(likely_systems, hit_system_tracker, sys_file):
    print(_outfile_header(), file=sys_file)
    if likely_systems:
        print("# Systems found:\n", file=sys_file)
        for system in likely_systems:
            sys_serializer = TxtLikelySystemSerializer()
            print(sys_serializer.serialize(system, hit_system_tracker), file=sys_file)
    else:
        print("# No Likely Systems found", file=sys_file)


def likely_systems_to_tsv(likely_systems, hit_system_tracker, sys_file):
    """
    print likely systems occurrences (from unordered replicon)
    in a file in human readable format

    :param likely_systems: list of systems found
    :type likely_systems: list of :class:`macsypy.system.LikelySystem` objects
    :param hit_system_tracker: a filled HitSystemTracker.
    :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :return: None
    """
    print(_outfile_header(), file=sys_file)
    if likely_systems:
        print("# Likely Systems found:\n", file=sys_file)
        print(TsvLikelySystemSerializer.header, file=sys_file)
        for l_system in likely_systems:
            sys_serializer = TsvLikelySystemSerializer()
            print(sys_serializer.serialize(l_system, hit_system_tracker), file=sys_file)
    else:
        print("# No Likely Systems found", file=sys_file)


def unlikely_systems_to_txt(unlikely_systems, sys_file):
    """
    print hits (from unordered replicon) which probably does not make a system occurrences
    in a file in human readable format

    :param unlikely_systems: list of :class:`macsypy.system.UnLikelySystem` objects
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :return: None
    """
    print(_outfile_header(), file=sys_file)
    if unlikely_systems:
        print("# Unlikely Systems found:\n", file=sys_file)
        for system in unlikely_systems:
            sys_serializer = TxtUnikelySystemSerializer()
            print(sys_serializer.serialize(system), file=sys_file)
            print("=" * 60, file=sys_file)
    else:
        print("# No Unlikely Systems found", file=sys_file)


def main(args=None, loglevel=None):
    """
    main entry point to MacSyFinder do some check before to launch :func:`main_search_systems` which is
    the real function that perform a search

    :param args: the arguments passed on the command line without the program name
    :type args: List of string
    :param loglevel: the output verbosity
    :type loglevel: a positive int or a string among 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    args = sys.argv[1:] if args is None else args
    parser, parsed_args = parse_args(args)

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
                raise ValueError(f"'{working_dir}' already exists and is not a empty")
        else:
            raise ValueError(f"'{working_dir}' already exists and is not a directory")

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
            parser.print_help()
            print()
            sys.tracebacklimit = 0
            raise OptionError("argument --models or --previous-run is required.")
        elif not parsed_args.previous_run and not parsed_args.sequence_db:
            parser.print_help()
            print()
            sys.tracebacklimit = 0
            raise OptionError("argument --sequence-db or --previous-run is required.")
        elif not parsed_args.previous_run and not parsed_args.db_type:
            parser.print_help()
            print()
            sys.tracebacklimit = 0
            raise OptionError("argument --db-type or --previous-run is required.")

        _log.info(f"command used: {' '.join(sys.argv)}")

        models = ModelBank()
        genes = GeneBank()
        profile_factory = ProfileFactory(config)

        logger.info("\n{:#^70}".format(" Searching systems "))
        all_systems, rejected_clusters = search_systems(config, models, genes, profile_factory, logger)

        track_multi_systems_hit = HitSystemTracker(all_systems)
        if config.db_type() in ('gembase', 'ordered_replicon'):
            #############################
            # Ordered/Gembase replicons #
            #############################

            ###########################
            # select the best systems #
            ###########################
            logger.info("\n{:#^70}".format(" Computing best solutions "))
            best_solutions = []
            one_best_solution = []

            # group systems found by replicon
            # before to search best system combination
            import time
            for rep_name, syst_group in itertools.groupby(all_systems, key=lambda s: s.replicon_name):
                syst_group = list(syst_group)
                logger.info(f"Computing best solutions for {rep_name} (nb of systems {len(syst_group)})")
                t0 = time.time()
                best_sol_4_1_replicon, score = find_best_solutions(syst_group)
                t1 = time.time()
                logger.info(f"It took {t1 - t0:.2f}sec to find best solution ({score:.2f}) for replicon {rep_name}")
                # if several solutions are equivalent same number of system and score is same
                # store all equivalent solution in best_solution => all_best_systems
                # pick one in one_best_solution => best_systems
                best_solutions.extend(best_sol_4_1_replicon)
                one_best_solution.append(best_sol_4_1_replicon[0])

            ##############################
            # Write the results in files #
            ##############################
            logger.info("\n{:#^70}".format(" Writing down results "))
            system_filename = os.path.join(config.working_dir(), "all_systems.txt")
            tsv_filename = os.path.join(config.working_dir(), "all_systems.tsv")

            with open(system_filename, "w") as sys_file:
                systems_to_txt(all_systems, track_multi_systems_hit, sys_file)

            with open(tsv_filename, "w") as tsv_file:
                systems_to_tsv(all_systems, track_multi_systems_hit, tsv_file)

            cluster_filename = os.path.join(config.working_dir(), "rejected_clusters.txt")
            with open(cluster_filename, "w") as clst_file:
                rejected_clusters.sort(key=lambda clst: (clst.replicon_name, clst.model, clst.hits))
                rejected_clst_to_txt(rejected_clusters, clst_file)
            if not (all_systems or rejected_clusters):
                logger.info("No Systems found in this dataset.")

            tsv_filename = os.path.join(config.working_dir(), "all_best_solutions.tsv")
            with open(tsv_filename, "w") as tsv_file:
                solutions_to_tsv(best_solutions, track_multi_systems_hit, tsv_file)

            tsv_filename = os.path.join(config.working_dir(), "best_solution.tsv")
            with open(tsv_filename, "w") as tsv_file:
                # flattern the list and sort it
                one_best_solution = [syst for sol in one_best_solution for syst in sol]
                one_best_solution.sort(key=lambda syst: (syst.replicon_name, syst.position[0], syst.model.fqn, - syst.score))
                systems_to_tsv(one_best_solution, track_multi_systems_hit, tsv_file)
        else:
            #######################
            # Unordered replicons #
            #######################

            ##############################
            # Write the results in files #
            ##############################
            logger.info("\n{:#^70}".format(" Writing down results "))

            system_filename = os.path.join(config.working_dir(), "all_systems.txt")
            with open(system_filename, "w") as sys_file:
                likely_systems_to_txt(all_systems, track_multi_systems_hit, sys_file)

            # forbidden = [s for s in all_systems if s.forbidden_occ]
            # system_filename = os.path.join(config.working_dir(), "forbidden_components.tsv")
            # with open(system_filename, "w") as sys_file:
            #     likely_systems_to_tsv(forbidden, track_multi_systems_hit, sys_file)

            system_filename = os.path.join(config.working_dir(), "all_systems.tsv")
            with open(system_filename, "w") as sys_file:
                likely_systems_to_tsv(all_systems, track_multi_systems_hit, sys_file)

            cluster_filename = os.path.join(config.working_dir(), "uncomplete_systems.txt")
            with open(cluster_filename, "w") as clst_file:
                unlikely_systems_to_txt(rejected_clusters, clst_file)

            if not (all_systems or rejected_clusters):
                logger.info("No Systems found in this dataset.")

    logger.info("END")


if __name__ == "__main__":
    main()
