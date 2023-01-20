#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2023  Institut Pasteur (Paris) and CNRS.           #
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

"""
Main entrypoint to macsyfinder
"""

import sys
import os
import argparse
import logging
import itertools
import signal
import time

from operator import attrgetter  # To be used with "sorted"
from textwrap import dedent

import colorlog

_log = colorlog.getLogger('macsypy')
import pandas as pd

import macsypy
from macsypy.config import MacsyDefaults, Config
from macsypy.cluster import Cluster
from macsypy.registries import ModelRegistry, scan_models_dir
from macsypy.definition_parser import DefinitionParser
from macsypy.search_genes import search_genes
from macsypy.database import Indexes, RepliconDB
from macsypy.error import OptionError, Timeout
from macsypy import cluster
from macsypy.hit import get_best_hits, HitWeight, MultiSystem, LonerMultiSystem, \
    sort_model_hits, compute_best_MSHit
from macsypy.system import OrderedMatchMaker, UnorderedMatchMaker, System, LikelySystem, UnlikelySystem, HitSystemTracker
from macsypy.utils import get_def_to_detect, get_replicon_names, parse_time
from macsypy.profile import ProfileFactory
from macsypy.model import ModelBank
from macsypy.gene import GeneBank
from macsypy.solution import find_best_solutions, combine_clusters, combine_multisystems
from macsypy.serialization import TxtSystemSerializer, TxtLikelySystemSerializer, TxtUnikelySystemSerializer, \
    TsvSystemSerializer, TsvSolutionSerializer, TsvLikelySystemSerializer, TsvSpecialHitSerializer, TsvRejectedCandidatesSerializer


def alarm_handler(signum, frame):
    _log.critical("Timeout is over. Aborting")
    for h in _log.handlers:
        h.flush()
    # I exit wit 0 otherwise in parallel_msf the job will be retry
    # on an other machine. we don't want that.
    #sys.exit(0)
    raise Timeout()


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
- Pandas {pd.__version__}

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
    model_dirs = config.models_dir()
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
    :rtype: :class:`argparse.Namespace` object.
    """
    parser = argparse.ArgumentParser(
        epilog="For more details, visit the MacSyFinder website and see the MacSyFinder documentation.",
        # formatter_class=ArgumentDefaultsHelpRawTextFormatter,
        formatter_class=argparse.RawTextHelpFormatter,
        description=dedent(r'''



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
(default: {msf_def['replicon_topology']})
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
    cut_ga_group = hmmer_options.add_mutually_exclusive_group()
    cut_ga_group.add_argument('--no-cut-ga',
                               action='store_true',
                               default=None,
                               help=f"""By default the MSF try to applied a threshold per profile by using the
hmmer -cut-ga option. This is possible only if the GA bit score is present in the profile otherwise 
MF switch to use the --e-value-search (-E in hmmsearch). 
If this option is set the --e-value-search option is used for all profiles regardless the presence of 
the a GA bit score in the profiles.
(default: {not msf_def['cut_ga']})""")
    cut_ga_group.add_argument('--cut-ga',
                               action='store_true',
                               default=None,
                               help=f"""By default the MSF try to applied a threshold per profile by using the
hmmer -cut-ga option. This is possible only if the GA bit score is present in the profile otherwise 
MSF switch to use the --e-value-search (-E in hmmsearch). 
But the modeler can override this default behavior to do not use cut_ga but --e-value-search instead (-E in hmmsearch).
The user can reestablish the general MSF behavior, be sure the profiles contain the GA bit score.
(default: {msf_def['cut_ga']})""")

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
                               help=f"""the weight of a accessory component in cluster scoring
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
    score_options.add_argument('--out-of-cluster',
                               action='store',
                               type=float,
                               default=None,
                               help=f"""the weight modifier for a hit which is a
 - true loner (not in cluster)
 - or multi-system (from an other system) 
(default:{msf_def['out_of_cluster_weight']})""")

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
    dir_options.add_argument('--index-dir',
                             action='store',
                             default=None,
                             help="Specifies the path to a directory to store/read the sequence index when the sequence-db dir is not writable.")
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
0 mean that all threads available will be used.
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

    general_options.add_argument("--timeout",
                                 action='store',
                                 default=None,
                                 type=parse_time,
                                 help="""In some case msf can take a long time to find the best solution (in 'gembase' and 'ordered_replicon mode').
The timeout is per replicon. If this step reach the timeout, the replicon is skipped (for gembase mode the analyse of other replicons continue).
NUMBER[SUFFIX]  NUMBER seconds. SUFFIX may be 's' for seconds (the default), 'm' for minutes, 'h' for hours or 'd' for days
for instance 1h2m3s means 1 hour 2 min 3 sec. NUMBER must be an integer.
""")

    parsed_args = parser.parse_args(args)
    if parsed_args.cfg_file and parsed_args.previous_run:
        # argparse does not allow to have mutually exclusive option  in a argument group
        # I prefer to have these 2 options in general options group
        # so I mimic the exclusive_group behavior
        parser.print_usage()
        print("macsyfinder: error: argument --previous-run: not allowed with argument --cfg-file")
        sys.exit(2)

    return parser, parsed_args


def search_systems(config, model_registry, models_def_to_detect, logger):
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
    :param model_registry: the registry of all models
    :type model_registry: :class:`macsypy.registries.ModelRegistry` object
    :param models_def_to_detect: the defintions to detect
    :type models_def_to_detect: list of :class:`macsypy.registries.DefinitionLocation` objects
    :param logger: The logger use to display information to the user.
                   It must be initialized. see :func:`macsypy.init_logger`
    :type logger: :class:`colorlog.Logger` object
    :return: the systems and rejected clusters found
    :rtype: ([:class:`macsypy.system.System`, ...], [:class:`macsypy.cluster.RejectedCAndidate`, ...])
    """
    working_dir = config.working_dir()
    config.save(path_or_buf=os.path.join(working_dir, config.cfg_name))

    # build indexes
    idx = Indexes(config)
    idx.build(force=config.idx())

    # create models
    model_bank = ModelBank()
    gene_bank = GeneBank()
    profile_factory = ProfileFactory(config)

    parser = DefinitionParser(config, model_bank, gene_bank, model_registry, profile_factory)
    parser.parse(models_def_to_detect)

    logger.info(f"MacSyFinder's results will be stored in working_dir{working_dir}")
    logger.info(f"Analysis launched on {config.sequence_db()} for model(s):")

    for model in models_def_to_detect:
        logger.info(f"\t- {model.fqn}")

    models_to_detect = [model_bank[model_loc.fqn] for model_loc in models_def_to_detect]
    all_genes = []
    for model in models_to_detect:
        genes = model.mandatory_genes + model.accessory_genes + model.neutral_genes + model.forbidden_genes
        # Exchangeable (formerly homologs/analogs) are also added because they can "replace" an important gene...
        ex_genes = []
        for m_gene in genes:
            ex_genes += m_gene.exchangeables
        all_genes += (genes + ex_genes)
    #############################################
    # this part of code is executed in parallel
    #############################################
    try:
        all_reports = search_genes(all_genes, config)
    except Exception as err:
        raise err
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
            systems, rejected_candidates = _search_in_ordered_replicon(hits_by_replicon, models_to_detect,
                                                                                       config, logger)
            return systems, rejected_candidates
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
    """

    :param hits_by_replicon:
    :param models_to_detect:
    :param config:
    :param logger:
    :return:
    """
    all_systems = []
    all_rejected_candidates = []
    rep_db = RepliconDB(config)
    
    for rep_name in hits_by_replicon:
        logger.info(f"\n{f' Hits analysis for replicon {rep_name} ':#^60}")
        rep_info = rep_db[rep_name]
        for model in models_to_detect:
            one_model_systems = []
            one_model_rejected_candidates = []
            logger.info(f"Check model {model.fqn}")
            # model.filter filter hit but also cast them in ModelHit
            mhits_related_one_model = model.filter(hits_by_replicon[rep_name])
            logger.debug(f"{f' hits related to {model.name} ':#^80}")
            hit_header_str = "id\trep_name\tpos\tseq_len\tgene_name\ti_eval\tscore\tprofile_cov\tseq_cov\tbeg_match\tend_match"
            hits_str = "".join([str(h) for h in mhits_related_one_model])
            logger.debug(f"\n{hit_header_str}\n{hits_str}")
            logger.debug("#" * 80)
            logger.info("Building clusters")
            hit_weights = HitWeight(**config.hit_weights())
            true_clusters, true_loners = cluster.build_clusters(mhits_related_one_model, rep_info, model, hit_weights)
            logger.debug(f"{' CLUSTERS ':#^80}")
            logger.debug("\n" + "\n".join([str(c) for c in true_clusters]))
            logger.debug(f"{' LONERS ':=^50}")
            logger.debug("\n" + "\n".join([str(c) for c in true_loners.values() if c.loner]))
            # logger.debug("{:=^50}".format(" MULTI-SYSTEMS hits "))
            # logger.debug("\n" + "\n".join([str(c.hits[0]) for c in special_clusters.values() if c.multi_system]))
            logger.debug("#" * 80)
            logger.info("Searching systems")
            clusters_combination = combine_clusters(true_clusters, true_loners, multi_loci=model.multi_loci)
            for one_clust_combination in clusters_combination:
                ordered_matcher = OrderedMatchMaker(model, redundancy_penalty=config.redundancy_penalty())
                res = ordered_matcher.match(one_clust_combination)
                if isinstance(res, System):
                    one_model_systems.append(res)
                else:
                    one_model_rejected_candidates.append(res)

            ###############################
            # MultiSystem Hits Management #
            ###############################
            # get multi systems from existing systems #
            hit_encondig_multisystems = set()  # for the same model (in the loop)
            for one_sys in one_model_systems:
                hit_encondig_multisystems.update(one_sys.get_hits_encoding_multisystem())

            logger.debug(f"{' MultiSystems ':#^80}")
            logger.debug("\n" + "\n".join([str(c) for c in true_clusters]))
            # Cast these hits in MultiSystem/LonerMultiSystem
            multi_systems_hits = []
            for hit in hit_encondig_multisystems:
                if not hit.loner:
                    multi_systems_hits.append(MultiSystem(hit))
                else:
                    multi_systems_hits.append(LonerMultiSystem(hit))
            # choose the best one
            ms_per_function = sort_model_hits(multi_systems_hits)
            best_ms = compute_best_MSHit(ms_per_function)
            # check if among rejected clusters with the MS, they can be created a new system
            best_ms = [Cluster([ms], model, hit_weights) for ms in best_ms]
            new_clst_combination = combine_multisystems(one_model_rejected_candidates, best_ms)
            for one_clust_combination in new_clst_combination:
                ordered_matcher = OrderedMatchMaker(model, redundancy_penalty=config.redundancy_penalty())
                res = ordered_matcher.match(one_clust_combination)
                if isinstance(res, System):
                    one_model_systems.append(res)
                else:
                    one_model_rejected_candidates.append(res)
            all_systems.extend(one_model_systems)
            all_rejected_candidates.extend(one_model_rejected_candidates)
    if all_systems:
        all_systems.sort(key=lambda syst: (syst.replicon_name, syst.position[0], syst.model.fqn, - syst.score))

    if not rep_db.guess_if_really_gembase():
        _log.warning(
            f"Most of replicons contains only ONE sequence are you sure that '{config.sequence_db()}' is a 'gembase'.")
    return all_systems, all_rejected_candidates


def _search_in_unordered_replicon(hits_by_replicon, models_to_detect, logger):
    """

    :param hits_by_replicon:
    :param models_to_detect:
    :param logger:
    :return:
    """
    likely_systems = []
    rejected_hits = []
    for rep_name in hits_by_replicon:
        logger.info(f"\n{f' Hits analysis for replicon {rep_name} ':#^60}")
        for model in models_to_detect:
            logger.info(f"Check model {model.fqn}")
            hits_related_one_model = model.filter(hits_by_replicon[rep_name])
            logger.debug("{:#^80}".format(" hits related to {} \n".format(model.name)))
            logger.debug("id\trep_name\tpos\tseq_len\tgene_name\ti_eval\tscore\tprofile_cov\tseq_cov\tbeg_match\tend_match")
            logger.debug("".join([str(h) for h in hits_related_one_model]))
            logger.debug("#" * 80)
            logger.info("Searching systems")
            hits_related_one_model = model.filter(hits_by_replicon[rep_name])
            if hits_related_one_model:
                unordered_matcher = UnorderedMatchMaker(model)
                res = unordered_matcher.match(hits_related_one_model)
                if isinstance(res, LikelySystem):
                    likely_systems.append(res)
                elif isinstance(res, UnlikelySystem):
                    rejected_hits.append(res)
                else:
                    logger.info(f"No hits related to {model.fqn } found.")
            else:
                logger.info(f"No hits found for model {model.fqn}")
    if likely_systems:
        likely_systems.sort(key=lambda syst: (syst.replicon_name, syst.position[0], syst.model.fqn))
    return likely_systems, rejected_hits


def _outfile_header(models_fam_name, models_version, skipped_replicons=None):
    """
    :return: The 2 first lines of each result file
    :rtype: str
    """
    header = f"""# macsyfinder {macsypy.__version__}
# models : {models_fam_name}-{models_version}
# {' '.join(sys.argv)}"""
    if skipped_replicons:
        header += "\n#"
        for rep_name in skipped_replicons:
            header += f"\n# WARNING: The replicon '{rep_name}' has been SKIPPED. Cannot be solved before timeout."
        header += "\n#"
    return header


def systems_to_tsv(models_fam_name, models_version, systems, hit_system_tracker, sys_file, skipped_replicons=None):
    """
    print systems occurrences in a file in tabulated  format

    :param systems: list of systems found
    :type systems: list of :class:`macsypy.system.System` objects
    :param hit_system_tracker: a filled HitSystemTracker.
    :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :param skipped_replicons: the replicons name for which msf reach the timeout
    :type skipped_replicons: list of str
    :return: None
    """
    print(_outfile_header(models_fam_name, models_version, skipped_replicons=skipped_replicons), file=sys_file)
    if systems:
        print("# Systems found:", file=sys_file)
        print(TsvSystemSerializer.header, file=sys_file)
        for system in systems:
            sys_serializer = TsvSystemSerializer()
            print(sys_serializer.serialize(system, hit_system_tracker), file=sys_file)
        warnings = _loner_warning(systems)
        if warnings:
            print("\n".join(warnings), file=sys_file)
    else:
        print("# No Systems found", file=sys_file)


def systems_to_txt(models_fam_name, models_version, systems, hit_system_tracker, sys_file, skipped_replicons=None):
    """
    print systems occurrences in a file in human readable format

    :param systems: list of systems found
    :type systems: list of :class:`macsypy.system.System` objects
    :param hit_system_tracker: a filled HitSystemTracker.
    :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :param skipped_replicons: the replicons name for which msf reach the timeout
    :type skipped_replicons: list of str
    :return: None
    """

    print(_outfile_header(models_fam_name, models_version, skipped_replicons=skipped_replicons), file=sys_file)
    if systems:
        print("# Systems found:\n", file=sys_file)
        for system in systems:
            sys_serializer = TxtSystemSerializer()
            print(sys_serializer.serialize(system, hit_system_tracker), file=sys_file)
            print("=" * 60, file=sys_file)

        warnings = _loner_warning(systems)
        if warnings:
            print("\n".join(warnings), file=sys_file)
    else:
        print("# No Systems found", file=sys_file)


def solutions_to_tsv(models_fam_name, models_version, solutions, hit_system_tracker, sys_file, skipped_replicons=None):
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
    :param skipped_replicons: the replicons name for which msf reach the timeout
    :type skipped_replicons: list of str
    :return: None
    """
    print(_outfile_header(models_fam_name, models_version, skipped_replicons=skipped_replicons), file=sys_file)
    if solutions:
        sol_serializer = TsvSolutionSerializer()
        print("# Systems found:", file=sys_file)
        print(sol_serializer.header, file=sys_file)

        for sol_id, solution in enumerate(solutions, 1):
            print(sol_serializer.serialize(solution, sol_id, hit_system_tracker), file=sys_file, end='')

            warnings = _loner_warning(solution.systems)
            if warnings:
                print("\n".join(warnings) + "\n", file=sys_file)
    else:
        print("# No Systems found", file=sys_file)


def _loner_warning(systems):
    """
    :param systems: sequence of systems
    :return: warning for loner which have less occurrences than systems occurrences in which this lone is used
             except if the loner is also multi system
    :rtype: list of string
    """
    warnings = []
    loner_tracker = {}
    for syst in systems:
        loners = syst.get_loners()
        for loner in loners:
            if loner.multi_system:
                # the loner multi_systems can appear in several systems
                continue
            elif loner in loner_tracker:
                loner_tracker[loner].append(syst)
            else:
                loner_tracker[loner] = [syst]
    for loner, systs in loner_tracker.items():
        if len(loner) < len(systs):
            # len(loners) count the number of loner occurrence the loner and its counterpart
            warnings.append(f"# WARNING Loner: there is only {len(loner)} occurrence(s) of loner '{loner.gene.name}' "
                            f"and {len(systs)} potential systems [{', '.join([s.id for s in systs])}]")

    return warnings


def summary_best_solution(models_fam_name, models_version, best_solution_path, sys_file, models_fqn, replicon_names,
                          skipped_replicons=None):
    """
    do a summary of best_solution in best_solution_path and write it on out_path
    a summary compute the number of system occurrence for each model and each replicon
    .. code-block:: text

        replicon        model_fqn_1  model_fqn_2  ....
        rep_name_1           1           2
        rep_name_2           2           0

    columns are separated by \t character

    :param str best_solution_path: the path to the best_solution file in tsv format
    :param sys_file: the file where to save the summary
    :param models_fqn: the fully qualified names of the models
    :type models_fqn: list of string
    :param replicon_names: the name of the replicons used
    :type replicon_names: list of string
    :param skipped_replicons: the replicons name for which msf reach the timeout
    :type skipped_replicons: list of str
    """
    skipped_replicons = skipped_replicons if skipped_replicons else set()
    print(_outfile_header(models_fam_name, models_version, skipped_replicons=skipped_replicons), file=sys_file)

    def fill_replicon(summary):
        """
        add row with 0 for all models for lacking replicons

        :param summary: the
        :type summary: :class:`pandas.DataFrame` object
        :return:
        :rtype: :class:`pandas.DataFrame` object
        """
        index_name = summary.index.name
        computed_replicons = set(summary.index)
        lacking_replicons = set(replicon_names) - computed_replicons - set(skipped_replicons)
        lacking_replicons = sorted(lacking_replicons)
        rows = pd.DataFrame({models: [0 * len(lacking_replicons)] for models in summary.columns}, index=lacking_replicons)
        summary = pd.concat([summary, rows], ignore_index=False)
        summary.index.name = index_name
        return summary

    def fill_models(summary):
        """
        add columns for lacking models (it means no occurence found)

        :param summary:
        :type summary: :class:`pandas.DataFrame` object
        :return:
        :rtype: :class:`pandas.DataFrame` object
        """
        computed_models = set(summary.columns)
        lacking_models = set(models_fqn) - computed_models
        lacking_models = sorted(lacking_models)
        for model in lacking_models:
            summary[model] = [0 for _ in summary.index]
        return summary

    try:
        best_solution = pd.read_csv(best_solution_path, sep='\t', comment='#')
    except pd.errors.EmptyDataError:
        # No results Found
        # may be there is no results so I have to report
        # may be the solution cannot be found So I do not to report (Warning)
        # may be the both one replicon without results one replicon not solved
        # So I have to report only the replicon without results
        replicon_to_report = list(set(replicon_names) - set(skipped_replicons))
        summary = pd.DataFrame(0, index=replicon_to_report, columns=models_fqn)
        summary.index.name = 'replicon'
    else:
        selection = best_solution[['replicon', 'sys_id', 'model_fqn']]
        dropped = selection.drop_duplicates(subset=['replicon', 'sys_id'])
        summary = pd.crosstab(index=dropped.replicon, columns=dropped['model_fqn'])
        summary = fill_replicon(summary)
        summary = fill_models(summary)

    summary.to_csv(sys_file, sep='\t')


def loners_to_tsv(models_fam_name, models_version, systems, sys_file):
    """
    get loners from valid systems and save them on file

    :param systems: the systems from which the loners are extract
    :type systems: list of :class:`macsypy.system.System` object
    :param sys_file: the file where loners are saved
    :type sys_file: file object open in write mode
    """
    print(_outfile_header(models_fam_name, models_version), file=sys_file)
    if systems:
        best_loners = set()
        for syst in systems:
            best_loners.update(syst.get_loners())
        if best_loners:
            serializer = TsvSpecialHitSerializer()
            loners = serializer.serialize(best_loners)
            print("# Loners found:", file=sys_file)
            print(loners, file=sys_file)
        else:
            print("# No Loners found", file=sys_file)
    else:
        print("# No Loners found", file=sys_file)


def multisystems_to_tsv(models_fam_name, models_version, systems, sys_file):
    """
    get multisystems from valid systems and save them on file

    :param systems: the systems from which the loners are extract
    :type systems: list of :class:`macsypy.system.System` object
    :param sys_file: the file where multisystems are saved
    :type sys_file: file object open in write mode
    """
    print(_outfile_header(models_fam_name, models_version), file=sys_file)
    if systems:
        best_multisystems = set()
        for syst in systems:
            best_multisystems.update(syst.get_multisystems())
        if best_multisystems:
            serializer = TsvSpecialHitSerializer()
            multisystems = serializer.serialize(best_multisystems)
            print("# Multisystems found:", file=sys_file)
            print(multisystems, file=sys_file)
        else:
            print("# No Multisystems found", file=sys_file)
    else:
        print("# No Multisystems found", file=sys_file)


def rejected_candidates_to_txt(models_fam_name, models_version, rejected_candidates, cand_file,
                               skipped_replicons=None):
    """
    print rejected clusters in a file

    :param rejected_candidates: list of candidates which does not contitute a system
    :type rejected_candidates: list of :class:`macsypy.system.RejectedCandidate` objects
    :param cand_file: The file where to write down the rejected candidates
    :type cand_file: file object
    :param skipped_replicons: the replicons name for which msf reach the timeout
    :type skipped_replicons: list of str
    :return: None
    """
    print(_outfile_header(models_fam_name, models_version, skipped_replicons=skipped_replicons), file=cand_file)
    if rejected_candidates:
        print("# Rejected candidates:\n", file=cand_file)
        for rej_cand in rejected_candidates:
            print(rej_cand, file=cand_file, end='')
            print("=" * 60, file=cand_file)
    else:
        print("# No Rejected candidates", file=cand_file)


def rejected_candidates_to_tsv(models_fam_name, models_version, rejected_candidates, cand_file,
                               skipped_replicons=None):
    """
    print rejected clusters in a file

    :param rejected_candidates: list of candidates which does not contitute a system
    :type rejected_candidates: list of :class:`macsypy.system.RejectedCandidate` objects
    :param cand_file: The file where to write down the rejected candidates
    :type cand_file: file object
    :param skipped_replicons: the replicons name for which msf reach the timeout
    :type skipped_replicons: list of str
    :return: None
    """
    print(_outfile_header(models_fam_name, models_version, skipped_replicons=skipped_replicons), file=cand_file)
    if rejected_candidates:
        serializer = TsvRejectedCandidatesSerializer()
        rej_candidates = serializer.serialize(rejected_candidates)
        print("# Rejected candidates found:", file=cand_file)
        print(rej_candidates, file=cand_file, end='')
    else:
        print("# No Rejected candidates", file=cand_file)


def likely_systems_to_txt(models_fam_name, models_version, likely_systems, hit_system_tracker, sys_file):
    """
    print likely systems occurrences (from unordered replicon)
    in a file in text human readable format
    :param likely_systems: list of systems found
    :type likely_systems: list of :class:`macsypy.system.LikelySystem` objects
    :param hit_system_tracker: a filled HitSystemTracker.
    :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
    :param sys_file: file object
    :return: None
    """
    print(_outfile_header(models_fam_name, models_version), file=sys_file)
    if likely_systems:
        print("# Systems found:\n", file=sys_file)
        for system in likely_systems:
            sys_serializer = TxtLikelySystemSerializer()
            print(sys_serializer.serialize(system, hit_system_tracker), file=sys_file)
    else:
        print("# No Likely Systems found", file=sys_file)


def likely_systems_to_tsv(models_fam_name, models_version, likely_systems, hit_system_tracker, sys_file):
    """
    print likely systems occurrences (from unordered replicon)
    in a file in tabulated separeted value (tsv) format

    :param likely_systems: list of systems found
    :type likely_systems: list of :class:`macsypy.system.LikelySystem` objects
    :param hit_system_tracker: a filled HitSystemTracker.
    :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :return: None
    """
    print(_outfile_header(models_fam_name, models_version), file=sys_file)
    if likely_systems:
        print("# Likely Systems found:\n", file=sys_file)
        print(TsvLikelySystemSerializer.header, file=sys_file)
        for l_system in likely_systems:
            sys_serializer = TsvLikelySystemSerializer()
            print(sys_serializer.serialize(l_system, hit_system_tracker), file=sys_file)
    else:
        print("# No Likely Systems found", file=sys_file)


def unlikely_systems_to_txt(models_fam_name, models_version, unlikely_systems, sys_file):
    """
    print hits (from unordered replicon) which probably does not make a system occurrences
    in a file in human readable format

    :param unlikely_systems: list of :class:`macsypy.system.UnLikelySystem` objects
    :param sys_file: The file where to write down the systems occurrences
    :type sys_file: file object
    :return: None
    """
    print(_outfile_header(models_fam_name, models_version), file=sys_file)
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

    if parsed_args.list_models:
        print(list_models(parsed_args), file=sys.stdout)
        sys.exit(0)

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

    #############################
    # command seems Ok Let's go #
    #############################
    _log.info(get_version_message())
    _log.info(f"command used: {' '.join(sys.argv)}")

    ########################################
    # compute which model I have to search #
    ########################################
    model_registry = ModelRegistry()
    for model_dir in config.models_dir():
        try:
            models_loc_available = scan_models_dir(model_dir,
                                                   profile_suffix=config.profile_suffix(),
                                                   relative_path=config.relative_path())
            for model_loc in models_loc_available:
                model_registry.add(model_loc)
        except PermissionError as err:
            _log.warning(f"{model_dir} is not readable: {err} : skip it.")

    try:
        models_def_to_detect, models_fam_name, models_version = get_def_to_detect(config.models(), model_registry)
    except KeyError as err:
        sys.exit(f"macsyfinder: {err}")
    _log.info(f"\nmodels used: {models_fam_name}-{models_version}")
    logger.info(f"\n{f' Searching systems ':#^70}")
    all_systems, rejected_candidates = search_systems(config, model_registry, models_def_to_detect, logger)
    track_multi_systems_hit = HitSystemTracker(all_systems)
    skipped_replicons = []
    
    if config.db_type() in ('gembase', 'ordered_replicon'):
        #############################
        # Ordered/Gembase replicons #
        #############################

        ###########################
        # select the best systems #
        ###########################
        logger.info(f"\n{f' Computing best solutions ':#^70}")
        all_best_solutions = []
        one_best_solution = []

        # group systems found by replicon
        # before to search best system combination
        for rep_name, syst_group in itertools.groupby(all_systems, key=lambda s: s.replicon_name):
            syst_group = list(syst_group)
            logger.info(f"Computing best solutions for {rep_name} (nb of candidate systems {len(syst_group)})")

            timeout = config.timeout()
            if timeout:
                # in some case best_solution take too much time
                # user can define a timeout by default set to 0
                signal.signal(signal.SIGALRM, alarm_handler)
                signal.alarm(config.timeout())
                _log.debug(f"set time out to {timeout} sec.")
            try:
                find_best_solutions_start = time.perf_counter()
                best_sol_4_1_replicon, score = find_best_solutions(syst_group)
                find_best_solutions_stop = time.perf_counter()
            except Timeout:
                _log.error(f"The {rep_name} cannot be solved in time skip it!")
                skipped_replicons.append(rep_name)
                continue
            if timeout:
                _log.debug("Cancel the time out.")
                signal.alarm(0)
                signal.signal(signal.SIGALRM, signal.SIG_DFL)

            logger.info(f"It took {find_best_solutions_stop - find_best_solutions_start:.2f}sec to find best solution"
                        f" ({score:.2f}) for replicon {rep_name}")
            # if several solutions are equivalent same number of system and score is same
            # store all equivalent solution in all_best_solution => all_best_systems
            # pick one in one_best_solution => best_systems
            all_best_solutions.extend(best_sol_4_1_replicon)
            one_best_solution.append(best_sol_4_1_replicon[0])

        ##############################
        # Write the results in files #
        ##############################
        logger.info(f"""\n{f" Writing down results in '{os.path.basename(config.working_dir())}' ":#^70}""")
        system_filename = os.path.join(config.working_dir(), "all_systems.txt")
        tsv_filename = os.path.join(config.working_dir(), "all_systems.tsv")

        with open(system_filename, "w") as sys_file:
            systems_to_txt(models_fam_name, models_version, all_systems, track_multi_systems_hit, sys_file,
                           skipped_replicons=skipped_replicons)

        with open(tsv_filename, "w") as tsv_file:
            systems_to_tsv(models_fam_name, models_version, all_systems, track_multi_systems_hit, tsv_file,
                           skipped_replicons=skipped_replicons)

        cluster_filename = os.path.join(config.working_dir(), "rejected_candidates.txt")
        with open(cluster_filename, "w") as clst_file:
            rejected_candidates.sort(key=lambda clst: (clst.replicon_name, clst.model, clst.hits))
            rejected_candidates_to_txt(models_fam_name, models_version, rejected_candidates, clst_file)
        if not (all_systems or rejected_candidates):
            logger.info("No Systems found in this dataset.")

        cluster_filename = os.path.join(config.working_dir(), "rejected_candidates.tsv")
        with open(cluster_filename, "w") as clst_file:
            rejected_candidates_to_tsv(models_fam_name, models_version, rejected_candidates, clst_file,
                                       skipped_replicons=skipped_replicons)

        tsv_filename = os.path.join(config.working_dir(), "all_best_solutions.tsv")
        with open(tsv_filename, "w") as tsv_file:
            solutions_to_tsv(models_fam_name, models_version, all_best_solutions, track_multi_systems_hit, tsv_file,
                             skipped_replicons=skipped_replicons)

        best_solution_filename = os.path.join(config.working_dir(), "best_solution.tsv")
        with open(best_solution_filename, "w") as best_solution_file:
            one_best_solution = [syst for sol in one_best_solution for syst in sol]
            one_best_solution.sort(key=lambda syst: (syst.replicon_name, syst.position[0], syst.model.fqn, - syst.score))
            systems_to_tsv(models_fam_name, models_version, one_best_solution, track_multi_systems_hit, best_solution_file,
                           skipped_replicons=skipped_replicons)

        loners_filename = os.path.join(config.working_dir(), "best_solution_loners.tsv")
        with open(loners_filename, "w") as loners_file:
            loners_to_tsv(models_fam_name, models_version, one_best_solution, loners_file)

        multisystems_filename = os.path.join(config.working_dir(), "best_solution_multisystems.tsv")
        with open(multisystems_filename, "w") as multisystems_file:
            multisystems_to_tsv(models_fam_name, models_version, one_best_solution, multisystems_file)

        summary_filename = os.path.join(config.working_dir(), "best_solution_summary.tsv")
        with open(summary_filename, "w") as summary_file:
            models_fqn = [m.fqn for m in models_def_to_detect]
            replicons_names = get_replicon_names(config.sequence_db(), config.db_type())
            summary_best_solution(models_fam_name, models_version, best_solution_filename, summary_file, models_fqn, replicons_names,
                                  skipped_replicons=skipped_replicons)

    else:
        #######################
        # Unordered replicons #
        #######################

        ##############################
        # Write the results in files #
        ##############################
        logger.info(f"""\n{f" Writing down results in '{os.path.basename(config.working_dir())}' ":#^70}""")

        system_filename = os.path.join(config.working_dir(), "all_systems.txt")
        with open(system_filename, "w") as sys_file:
            likely_systems_to_txt(models_fam_name, models_version, all_systems, track_multi_systems_hit, sys_file)

        # forbidden = [s for s in all_systems if s.forbidden_occ]
        # system_filename = os.path.join(config.working_dir(), "forbidden_components.tsv")
        # with open(system_filename, "w") as sys_file:
        #     likely_systems_to_tsv(forbidden, track_multi_systems_hit, sys_file)

        system_filename = os.path.join(config.working_dir(), "all_systems.tsv")
        with open(system_filename, "w") as sys_file:
            likely_systems_to_tsv(models_fam_name, models_version, all_systems, track_multi_systems_hit, sys_file)

        cluster_filename = os.path.join(config.working_dir(), "uncomplete_systems.txt")
        with open(cluster_filename, "w") as clst_file:
            unlikely_systems_to_txt(models_fam_name, models_version, rejected_candidates, clst_file)

        if not (all_systems or rejected_candidates):
            logger.info("No Systems found in this dataset.")
    if skipped_replicons:
        for rep_name in skipped_replicons:
            _log.error(f"The replicon {rep_name} cannot be solved before timeout. SKIP IT.")
    logger.info("END")


if __name__ == "__main__":
    main()
