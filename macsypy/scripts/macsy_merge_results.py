#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2020  Institut Pasteur (Paris) and CNRS.           #
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
import sys
import macsypy
import argparse
from typing import List

import colorlog

from macsypy.config import MacsyDefaults


def merge_files(files, out, ignore=None, keep_first=None, header=""):
    """

    :param files: the list of files to merge
    :type files: list of str
    :param str out: the path to the merged file
    :param str ignore: a string which start the lines to ignore
    :param str keep_first: a string which start the line which must be keep
                       only the first time
    :param str header: The header of the merged file
    :return:
    """
    with open(out, 'w') as f_out:
        f_out.write(header)
        if keep_first is None:
            first = False
        else:
            first = True
        for file in files:
            _log.debug(f"Merging {file}")
            with open(file) as fh_in:
                for line in fh_in:
                    if first and line.startswith(keep_first):
                        first = False
                    elif ignore and line.startswith(ignore):
                        continue
                    elif keep_first and line.startswith(keep_first):
                        continue

                    f_out.write(line)


def merge_and_reindex(files, out, header=""):
    """
    merge all_best_solutions and reindex the sol_id column

    :param files: the list of files to merge
    :type files: list of str
    :param str out: the path to the merged file
    """
    with open(out, 'w') as f_out:
        f_out.write(header)
        last_sol_id = 0
        first = True
        for file in files:
            _log.debug(f"Merging {file}")
            with open(file) as fh_in:
                for line in fh_in:
                    if first and line.startswith('sol_id'):
                        first = False
                        new_line = line
                    elif line.startswith('sol_id'):
                        continue
                    elif line.startswith('#'):
                        continue
                    elif line.startswith('\n'):
                        new_line = line
                    elif line:
                        fields = line.split('\t')
                        try:
                            new_sol_id = int(fields[0]) + last_sol_id
                        except ValueError as err:
                            _log.critical(f"Cannot reindex int({fields[0]}) + {last_sol_id}: {err}")
                            raise ValueError() from None
                        fields[0] = str(new_sol_id)
                        new_line = '\t'.join(fields)
                    f_out.write(new_line)
            last_sol_id = new_sol_id


def merge_results(results_dirs, out_dir='.'):
    """

    :param results_dirs: The list of macsyfinder results directories to merge
    :type results_dirs: list of str
    :param str out_dir: the path to the directory where to store the merged files
    """
    out_file = os.path.join(out_dir, 'merged_all_best_solutions.tsv')
    _log.info(f"Merging 'all_best_solutions.tsv' in to '{out_file}'")
    all_best_solutions_files = [os.path.join(d, 'all_best_solutions.tsv') for d in results_dirs]
    header = f"""# parallel_msf {macsypy.__version__}
# merged all_best_solutions.tsv
# systems found:
"""
    merge_and_reindex(all_best_solutions_files, out_file, header=header)

    out_file = os.path.join(out_dir, 'merged_best_solution.tsv')
    _log.info(f"Merging 'best_solution.tsv' in to '{out_file}'")
    best_solution_files = [os.path.join(d, 'best_solution.tsv') for d in results_dirs]
    header = f"""# parallel_msf {macsypy.__version__}
# merged best_solution.tsv
# systems found:
"""
    merge_files(best_solution_files, out_file, ignore="#", keep_first="replicon", header=header)

    out_file = os.path.join(out_dir, "merged_all_systems.tsv")
    _log.info(f"Merging 'all_systems.tsv' in to '{out_file}'")
    all_systems_files = [os.path.join(d, 'all_systems.tsv') for d in results_dirs]
    header = f"""# parallel_msf {macsypy.__version__}
# merged all_systems.tsv
# systems found:
"""
    merge_files(all_systems_files, out_file, ignore="#", keep_first="replicon", header=header)

    out_file = os.path.join(out_dir, 'merged_all_systems.txt')
    _log.info(f"Merging 'all_systems.txt' in to '{out_file}'")
    all_systems_files = [os.path.join(d, 'all_systems.txt') for d in results_dirs]
    header = f"""# parallel_msf {macsypy.__version__}
# merged all_systems.txt
# systems found:
"""
    merge_files(all_systems_files, out_file, ignore="#", header=header)

    out_file = os.path.join(out_dir, 'merged_rejected_clusters.txt')
    _log.info(f"Merging 'rejected_clusters.txt' in to '{out_file}'")
    all_rejected_files = [os.path.join(d, 'rejected_clusters.txt') for d in results_dirs]
    header = f"""# parallel_msf {macsypy.__version__}
# merged rejected_clusters.txt
# Rejected clusters:
"""
    merge_files(all_rejected_files, out_file, ignore="#", header=header)


def parse_args(args:  List[str]) -> argparse.Namespace:
    """

    :param args: the arguments passed on the command line without the first elemnet
    :type args: list of str
    :return: the command line and options arguments parsed
    :rtype: :class:`argparse.Namespace` object
    """
    description = """Merge the different files from several macsyfinder results in one.
    
    - merge the 'best_solution.tsv' in to 'merged_best_solution.tsv'
    - merge the 'all_best_solutions.tsv' in to `merged_all_best_solutions'
    - merge the 'all_systems.tsv' in to 'merged_all_systems.tsv'
    - merge the 'all_systems.txt' in to 'merged_all_systems.txt'
    - merge the 'rejected_clusters.txt' in to 'merged_rejected_clusters.txt'
"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=description)
    parser.add_argument('results_dirs',
                         nargs='+',
                         help='Path to the macsyfinder results directories to merge eg : path/to/macsyfinder-date-hour')
    parser.add_argument('-o', '--outdir',
                        default='.',
                        help='The path to the directory where to write merged files.')
    parser.add_argument("--mute",
                        action='store_true',
                        default=False,
                        help="mute the log on stdout."
                             "(continue to log on macsy_merge_results.out)")

    verbosity_grp = parser.add_argument_group()
    verbosity_grp.add_argument('-v', '--verbose',
                               action='count',
                               default=0,
                               help='Increase verbosity of output (can be cumulative : -vv)')
    verbosity_grp.add_argument('-q', '--quiet',
                               action='count',
                               default=0,
                               help='Decrease verbosity of output (can be cumulative : -qq)'
                               )
    parsed_args = parser.parse_args(args)
    return parsed_args


def main(args=None, log_level=None) -> None:
    """
    main entry point to macsy_merge_results

    :param args: the arguments passed on the command line
    :type args: list of str
    :param log_level: the output verbosity
    :type log_level: a positive int or a string among 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    global _log

    args = sys.argv[1:] if args is None else args
    parsed_args = parse_args(args)

    outdir_err = None
    if not os.path.exists(parsed_args.outdir):
        try:
            os.mkdir(parsed_args.outdir)
        except PermissionError as err:
            outdir_err = f"Cannot create {parsed_args.outdir} : {err}"
    elif not os.path.isdir(parsed_args.outdir):
        outdir_err = f"{parsed_args.outdir} is not a directory"
    elif not os.access(parsed_args.outdir, os.W_OK):
        outdir_err = f"{parsed_args.outdir} is not writable"
    if outdir_err:
        log = colorlog.getLogger('macsypy.merge')
        stdout_handler = colorlog.StreamHandler(sys.stdout)
        stdout_formatter = colorlog.ColoredFormatter("%(log_color)s%(message)s",
                                                     datefmt=None,
                                                     reset=True,
                                                     log_colors={'CRITICAL': 'bold_red'},
                                                     secondary_log_colors={},
                                                     style='%'
                                                     )
        stdout_handler.setFormatter(stdout_formatter)
        log.addHandler(stdout_handler)
        log.critical(outdir_err)
        sys.tracebacklimit = 0
        raise IOError() from None

    macsypy.init_logger(log_file=os.path.join(parsed_args.outdir, 'macsy_merge_results.out'),
                        out=not parsed_args.mute)
    _log = colorlog.getLogger('macsypy.merge')

    if not log_level:
        # logs are specify from args options
        config = MacsyDefaults()
        log_level = max(config.log_level - (10 * parsed_args.verbose) + (10 * parsed_args.quiet), 1)
        macsypy.logger_set_level(log_level)
    else:
        # used by unit tests to mute or unmute logs
        macsypy.logger_set_level(log_level)

    merge_results(parsed_args.results_dirs, out_dir=parsed_args.outdir)


if __name__ == '__main__':
    main()
