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
import glob
import argparse
from typing import List, Dict
from itertools import groupby
from dataclasses import dataclass
from textwrap import dedent
import logging

import colorlog

import macsypy
from macsypy.config import MacsyDefaults, Config
from macsypy.database import Indexes
from macsypy.hit import get_best_hits

# _log is set in main func
_log = None


def get_version_message() -> str:
    """
    :return: the long description of the macsyfinder version
    :rtype: str
    """
    version = macsypy.__version__
    vers_msg = f"""macsyprofile {version}
Python {sys.version}

MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
See the COPYING file for details.

If you use this software please cite:
{macsypy.__citation__}
and don't forget to cite models used:
macsydata cite <model>
"""
    return vers_msg


def get_profile_len(path: str) -> int:
    """
    Parse the HMM profile to extract the length and the presence of GA bit threshold

    :param str path: The path to the hmm profile used to produced the hmm search output to analyse
    :return: the length, presence of ga bit threshold
    :rtype: tuple(int length, bool ga_threshold)
    """
    with open(path) as f:
        for l in f:
            if l.startswith("LENG"):
                length = int(l.split()[1])
                break
    return length


def get_gene_name(path: str, suffix: str) -> str:
    """

    :param str path: The path to the hmm output to analyse
    :param str suffix: the suffix of the hmm output file
    :return: the name of the analysed gene
    :rtype: str
    """
    f = os.path.basename(path)
    gene_name = f.replace(suffix, '')
    return gene_name


@dataclass
class LightHit:
    """
    Handle hmm hits
    """

    gene_name: str
    id: str
    seq_length: int
    replicon_name: str
    position: int
    i_eval: float
    score: float
    profile_coverage: float
    sequence_coverage: float
    begin_match: int
    end_match: int


    def __str__(self) -> str:
        return f"{self.id}\t{self.replicon_name}\t{self.position:d}\t{self.seq_length:d}\t{self.gene_name}\t" \
               f"{self.i_eval:.3e}\t{self.score:.3f}\t{self.profile_coverage:.3f}\t" \
               f"{self.sequence_coverage:.3f}\t{self.begin_match:d}\t{self.end_match:d}"


class HmmProfile:

    def __init__(self, gene_name, gene_profile_lg, hmmer_output, cfg):
        """
        :param gene: the gene corresponding to the profile search reported here
        :type gene: :class:`macsypy.gene.CoreGene` object
        :param hmmer_output: The path to the raw Hmmer output file
        :type hmmer_output: string
        :param cfg: the configuration object
        :type cfg: :class:`macsypy.config.Config` object
        """
        self.gene_name = gene_name
        self._hmmer_raw_out = hmmer_output
        self.gene_profile_lg = gene_profile_lg
        self.cfg = cfg


    def parse(self) -> List[LightHit]:
        """
        parse a hmm output file and extract all hits and do some basic computation (coverage profile)

        :return: The list of extracted hits
        """
        all_hits = []
        idx = Indexes(self.cfg)
        macsyfinder_idx = idx.build()
        my_db = self._build_my_db(self._hmmer_raw_out)
        self._fill_my_db(macsyfinder_idx, my_db)

        with open(self._hmmer_raw_out, 'r') as hmm_out:
            i_evalue_sel = self.cfg.i_evalue_sel()
            coverage_threshold = self.cfg.coverage_profile()
            hmm_hits = (x[1] for x in groupby(hmm_out, self._hit_start))
            # drop summary
            next(hmm_hits)
            for hmm_hit in hmm_hits:
                hit_id = self._parse_hmm_header(hmm_hit)
                seq_lg, position_hit = my_db[hit_id]

                replicon_name = self._get_replicon_name(hit_id)

                body = next(hmm_hits)
                h = self._parse_hmm_body(hit_id, self.gene_profile_lg, seq_lg, coverage_threshold,
                                         replicon_name, position_hit, i_evalue_sel, body)
                all_hits += h
            hits = sorted(all_hits, key=lambda h: - h.score)
        return hits


    def _build_my_db(self, hmm_output: str) -> Dict:
        """
        Build the keys of a dictionary object to store sequence identifiers of hits.

        :param hmm_output: the path to the hmmsearch output to parse.
        :type hmm_output: string
        :return: a dictionary containing a key for each sequence id of the hits
        :rtype: dict
        """
        d = {}
        with open(hmm_output) as hmm_file:
            hits = (x[1] for x in groupby(hmm_file, self._hit_start) if x[0])
            for h in hits:
                d[self._parse_hmm_header(h)] = None
        return d


    def _fill_my_db(self, macsyfinder_idx: str, db: Dict) -> None:
        """
        Fill the dictionary with information on the matched sequences

        :param macsyfinder_idx: the path the macsyfinder index corresponding to the dataset
        :type  macsyfinder_idx: string
        :param db: the database containing all sequence id of the hits.
        :type db: dict
        """
        with open(macsyfinder_idx, 'r') as idx:
            for l in idx:
                seqid, length, rank = l.split(';')
                if seqid in db:
                    db[seqid] = (int(length), int(rank))


    def _get_replicon_name(self, hit_id: str) -> str:

        replicon_name = {'unordered': 'unordered',
                         'ordered_replicon': 'UserReplicon',
                         'gembase': "_".join(hit_id.split('_')[:-1])}
        return replicon_name[self.cfg.db_type()]


    def _hit_start(self, line: str) -> bool:
        """
        :param line: the line to parse
        :type line: string
        :return: True if it's the beginning of a new hit in Hmmer raw output files.
         False otherwise
        :rtype: boolean.
        """
        return line.startswith(">>")


    def _parse_hmm_header(self, h_grp) -> str:
        """
        :param h_grp: the sequence of string return by groupby function representing the header of a hit
        :type h_grp: sequence of string (<itertools._grouper object at 0x7ff9912e3b50>)
        :returns: the sequence identifier from a set of lines that corresponds to a single hit
        :rtype: string
        """
        for line in h_grp:
            hit_id = line.split()[1]
        return hit_id


    def _parse_hmm_body(self, hit_id, gene_profile_lg, seq_lg, coverage_threshold, replicon_name,
                        position_hit, i_evalue_sel, b_grp):
        """
        Parse the raw Hmmer output to extract the hits, and filter them with threshold criteria selected
        ("coverage_profile" and "i_evalue_select" command-line parameters)

        :param str hit_id: the sequence identifier
        :param int gene_profile_lg: the length of the profile matched
        :paramint  seq_lg: the length of the sequence
        :param float coverage_threshold: the minimal coverage of the profile to be reached in the Hmmer alignment
                                        for hit selection.
        :param str replicon_name: the identifier of the replicon
        :param int position_hit: the rank of the sequence matched in the input dataset file
        :param float i_evalue_sel: the maximal i-evalue (independent evalue) for hit selection
        :param b_grp: the Hmmer output lines to deal with (grouped by hit)
        :type b_grp: list of list of strings
        :returns: a sequence of hits
        :rtype: list of :class:`macsypy.report.Hit` objects

        """
        first_line = next(b_grp)
        if not first_line.startswith('   #    score'):
            return []
        else:
            hits = []
            for line in b_grp:
                if line[0] == '\n':
                    return hits
                elif line.startswith(" ---   ------ ----- --------"):
                    pass
                else:
                    fields = line.split()
                    try:
                        # fields[2] = score
                        # fields[5] = i_evalue
                        # fields[6] = hmmfrom
                        # fields[7] = hmm to
                        # fields[9] = alifrom
                        # fields[10] = ali to
                        if len(fields) > 1:
                            _, _, score, _, _, i_evalue, hmm_from, hmm_to, _, ali_from, ali_to, *_ = fields
                            score = float(score)
                            i_evalue = float(i_evalue)
                            hmm_from = int(hmm_from)
                            hmm_to = int(hmm_to)
                            ali_from = int(ali_from)
                            ali_to = int(ali_to)
                            if i_evalue <= i_evalue_sel:
                                _log.debug(f"{hit_id} i_evalue {i_evalue} <= {i_evalue_sel} i_evalue_sel")
                                cov_profile = (hmm_to - hmm_from + 1) / gene_profile_lg
                                begin = int(fields[9])
                                end = int(fields[10])
                                cov_gene = (end - begin + 1) / seq_lg  # To be added in Gene: sequence_length
                                if cov_profile >= coverage_threshold:
                                    hits.append(LightHit(self.gene_name, hit_id, seq_lg, replicon_name, position_hit,
                                                         i_evalue, score, cov_profile, cov_gene, ali_from, ali_to))
                                    _log.debug(f"{hit_id} cov_profile {cov_profile} >= {coverage_threshold} "
                                               f"coverage_threshold: add hit")

                                else:
                                    _log.debug(f"{hit_id} cov_profile {cov_profile} < {coverage_threshold} "
                                               f"coverage_threshold: skip hit")
                            else:
                                _log.debug(f"{hit_id} i_evalue {i_evalue} > {i_evalue_sel} i_evalue_sel : skip hit")
                    except ValueError as err:
                        msg = f"Invalid line to parse :{line}:{err}"
                        _log.debug(msg)
                        raise ValueError(msg)


def header(cmd: List[str]) -> str:
    """

    :param cmd: the command use dto launch this analyse
    :return: The header of the result file
    """
    header = f"""# macsyprofile {macsypy.__version__}
# macsyprofile {' '.join(cmd)}
hit_id\treplicon_name\tposition_hit\thit_sequence_length\tgene_name\ti_eval\tscore\tprofile_coverage\tsequence_coverage\tbegin\tend"""
    return header


def init_logger(level='INFO', out=True):
    """

    :param level: The logger threshold could be a positive int or string
                  among: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'
    :param out: if the log message must be displayed
    :return: logger
    :rtype: :class:`logging.Logger` instance
    """
    logger = colorlog.getLogger('macsyprofile')
    logging = colorlog.logging.logging
    if isinstance(level, str):
        level = getattr(logging, level)
    if out:
        stdout_handler = colorlog.StreamHandler(sys.stderr)
        if level <= logging.DEBUG:
            msg_formatter = "%(log_color)s%(levelname)-8s : %(module)s: L %(lineno)d :%(reset)s %(message)s"
        else:
            msg_formatter = "%(log_color)s%(message)s"
        stdout_formatter = colorlog.ColoredFormatter(msg_formatter,
                                                     datefmt=None,
                                                     reset=True,
                                                     log_colors={
                                                         'DEBUG': 'cyan',
                                                         'INFO': 'green',
                                                         'WARNING': 'yellow',
                                                         'ERROR': 'red',
                                                         'CRITICAL': 'bold_red',
                                                     },
                                                     secondary_log_colors={},
                                                     style='%'
                                                     )
        stdout_handler.setFormatter(stdout_formatter)
        logger.addHandler(stdout_handler)
    else:
        null_handler = logging.NullHandler()
        logger.addHandler(null_handler)
    logger.setLevel(level)
    return logger


def verbosity_to_log_level(verbosity: int) -> int:
    """
    transform the number of -v option in loglevel
    :param int verbosity: number of -v option on the command line
    :return: an int corresponding to a logging level
    """
    level = max((logging.INFO - (10 * verbosity), 1))
    return level


def parse_args(args:  List[str]) -> argparse.Namespace:
    """

    :param args: The arguments provided on the command line
    :type args: List of strings [without the program name]
    :return: The arguments parsed
    :rtype: :class:`aprgparse.Namespace` object.
    """
    msf_def = MacsyDefaults()
    parser = argparse.ArgumentParser(
        epilog="For more details, visit the MacSyFinder website and see the MacSyFinder documentation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=dedent('''

         *            *               *                   * *       * 
    *           *               *   *   *  *    **                *  
      **     *    *   *  *     *                    *               *
         __  __    *       ____ *      ____        **   __ _ _  *    
        |  \/  | __ _  ___/ ___| _   _|  _ \ _ __ ___  / _(_) | ___ 
        | |\/| |/ _` |/ __\___ \| | | | |_) | '__/ _ \| |_| | |/ _ \ 
        | |  | | (_| | (__ ___) | |_| |  __/| | | (_) |  _| | |  __/
        |_|  |_|\__,_|\___|____/ \__, |_|   |_|  \___/|_| |_|_|\___|
                *                |___/    *                   *
     *      *   * *     *   **         *   *  *           *
      *      *         *        *    *              *
                 *                           *  *           *     *


    MacSyProfile - MacSyFinder profile helper tool
    '''))

    parser.add_argument('previous_run',
                        action='store',
                        help='The path to a macsyfinder results directory.'
                        )
    parser.add_argument('--coverage-profile',
                        action='store',
                        default=-1.,
                        type=float,
                        help=f"""Minimal profile coverage required for the hit alignment  with the profile to allow
the hit selection for systems detection. (default no threshold)"""
                        )
    parser.add_argument('--i-evalue-sel',
                        action='store',
                        type=float,
                        default=1.0e9,
                        help=f"""Maximal independent e-value for Hmmer hits to be selected for systems detection.
(default: no selection based on i-evalue)""")
    parser.add_argument('--best-hits',
                        choices=['score', 'i_eval', 'profile_coverage'],
                        action='store',
                        default=None,
                        help="If several hits match the same replicon, same gene. "
                             "Select only the best one (based on best 'score' or 'i_evalue' or 'profile_coverage')")
    parser.add_argument('-p', '--pattern',
                        action='store',
                        default='*',
                        help="pattern to filter the hmm files to analyse."
                        )
    parser.add_argument('-o', '--out',
                        action='store',
                        default=None,
                        help="the path to a file to write results.")
    parser.add_argument('-f', '--force',
                        action='store_true',
                        default=False,
                        help='force to write output even the file already exists (overwrite it).')
    parser.add_argument('-V', "--version",
                        action="version",
                        version=get_version_message())
    parser.add_argument("-v", "--verbosity",
                        action="count",
                        default=0,
                        help="""Increases the verbosity level. There are 4 levels:
Error messages (default), Warning (-v), Info (-vv) and Debug.(-vvv)""")
    parser.add_argument("--mute",
                        action="store_true",
                        default=False,
                        help=f"""Mute the log on stdout.
(continue to log on macsyfinder.log)
(default: {msf_def['mute']})""")

    parsed_args = parser.parse_args(args)

    return parsed_args


def main(args=None, log_level=None) -> None:
    """
    main entry point to macsyprofile

    :param args: the arguments passed on the command line without the program name
    :type args: List of string
    :param log_level: the output verbosity
    :type log_level: a positive int or a string among 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    global _log
    args = sys.argv[1:] if args is None else args
    parsed_args = parse_args(args)

    if log_level is None:
        log_level = verbosity_to_log_level(parsed_args.verbosity)
    _log = init_logger(log_level, out=(not parsed_args.mute))

    if not os.path.exists(parsed_args.previous_run):
        _log.critical(f"{parsed_args.previous_run}: No such directory.")
        sys.tracebacklimit = 0
        raise FileNotFoundError() from None
    elif not os.path.isdir(parsed_args.previous_run):
        _log.critical(f"{parsed_args.previous_run} is not a directory.")
        sys.tracebacklimit = 0
        raise ValueError() from None

    defaults = MacsyDefaults(i_evalue_sel=1.0e9, coverage_profile=-1.0)
    cfg = Config(defaults, parsed_args)

    msf_run_path = cfg.previous_run()
    hmmer_results = os.path.join(msf_run_path, cfg.hmmer_dir())
    hmm_suffix = cfg.res_search_suffix()
    profile_suffix = cfg.profile_suffix()
    if parsed_args.out:
        profile_report_path = os.path.normpath(parsed_args.out)
        dirname = os.path.normpath(os.path.dirname(parsed_args.out))
        if not os.path.exists(dirname):
            _log.critical(f"The {dirname} directory is not writable")
            sys.tracebacklimit = 0
            raise ValueError() from None
    else:
        profile_report_path = os.path.join(cfg.previous_run(), 'hmm_coverage.tsv')

    if os.path.exists(profile_report_path) and not parsed_args.force:
        _log.critical(f"The file {profile_report_path} already exists. "
                      f"Remove it or specify a new output name --out or use --force option")
        sys.tracebacklimit = 0
        raise ValueError() from None

    hmmer_files = sorted(glob.glob(os.path.join(hmmer_results, f"{parsed_args.pattern}{hmm_suffix}")))
    try:
        profiles_dir = os.path.join(cfg.models_dir(), cfg.models()[0], 'profiles')
    except IndexError:
        _log.critical(f"Cannot find models in conf file {msf_run_path}. "
                      f"May be these results have been generated with an old version of macsyfinder.")
        sys.tracebacklimit = 0
        raise ValueError() from None

    _log.debug(f"hmmer_files: {hmmer_files}")
    all_hits = []
    with open(profile_report_path, 'w') as prof_out:
        print(header(args), file=prof_out)
        for hmmer_out_path in hmmer_files:
            _log.info(f"parsing {hmmer_out_path}")
            gene_name = get_gene_name(hmmer_out_path, hmm_suffix)
            profile_path = os.path.join(profiles_dir, f"{gene_name}{profile_suffix}")
            gene_profile_len = get_profile_len(profile_path)
            hmm = HmmProfile(gene_name, gene_profile_len, hmmer_out_path, cfg)
            hits = hmm.parse()
            all_hits += hits
        if len(all_hits) > 0:
            if parsed_args.best_hits:
                # It's important to keep this sorting to have in last all_hits version
                # the hits with the same replicon_name and position sorted by score
                # the best score in first
                hits_by_replicon = {}
                for hit in all_hits:
                    if hit.replicon_name in hits_by_replicon:
                        hits_by_replicon[hit.replicon_name].append(hit)
                    else:
                        hits_by_replicon[hit.replicon_name] = [hit]
                all_hits = []
                for rep_name in hits_by_replicon:
                    hits_by_replicon[rep_name] = get_best_hits(hits_by_replicon[rep_name], key=parsed_args.best_hits)
                    all_hits += sorted(hits_by_replicon[rep_name], key=lambda h: h.position)

            all_hits = sorted(all_hits, key=lambda h: (h.gene_name, h.replicon_name, h.position, h.score))
            _log.info(f"found {len(all_hits)} hits")
            for hit in all_hits:
                print(hit, file=prof_out)
            _log.info(f"result is in '{profile_report_path}'")
        else:
            _log.info(f"No hit found")


if __name__ == '__main__':
    main()
