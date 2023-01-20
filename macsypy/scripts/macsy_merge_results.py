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

import os
import sys
import macsypy
import argparse
from typing import List

import colorlog
import pandas as pd

from macsypy.config import MacsyDefaults


def get_warning(path):
    """

    :param path: the path of the result file to parse
    :type path: str
    :return: the list of warning in the header
    :rtype: list of str
    """
    warn_to_report = []
    with open(path) as file:
        try:
            line = next(file)
            while line.startswith('#'):
                if line.startswith('# WARNING: The replicon'):
                    warn_to_report.append(line)
                line = next(file)
        except StopIteration:
            pass
    return warn_to_report


def merge_files(files: List[str], out: str, header: str,
                ignore: str = None, keep_first: str = None, skip_until=None) -> None:
    """

    :param files: the list of files to merge
    :type files: list of str
    :param str out: the path to the merged file
    :param str ignore: a string which start the lines to ignore
    :param str keep_first: a string which start the line which must be keep
                       only the first time
    :param skip_until: skip all lines until the condition is True
    :type skip_until: a fonction which test the line
    :param str header: The header of the merged file
    :return:
    """
    def get_header(result: bool, warnings):
        res_or_not = header if result else f"No {header}"
        header_str = f"""# parallel_msf {macsypy.__version__}
# merged {os.path.basename(files[0])}
# {res_or_not} found:
"""
        if warnings:
            header_str += '# \n'
            for warn in warnings:
                header_str += warn
            header_str += '# \n'

        return header_str
    warnings = []
    with open(out, 'w') as f_out:
        if keep_first is None:
            first = False
        else:
            first = True
        results = False
        for file in files:
            _log.debug(f"Merging {file}")
            warnings += get_warning(file)
            skip = bool(skip_until)
            with open(file) as fh_in:
                for line in fh_in:
                    if skip:
                        if skip_until(line):
                            skip = False
                        else:
                            continue
                    if first and line.startswith(keep_first):
                        first = False
                        f_out.write(get_header(True, warnings))
                    elif ignore and line.startswith(ignore):
                        continue
                    elif keep_first and line.startswith(keep_first):
                        continue
                    f_out.write(line)
                    results = True
                if not first:
                    f_out.write('\n')
        if not results:
            f_out.write(get_header(False, warnings) + '\n')


def merge_and_reindex(files: List[str], out: str,  header: str,
                      comment: str = None,  skip_until=None) -> None:
    """
    merge all_best_solutions and reindex the sol_id column

    :param files: the list of files to merge
    :type files: list of str
    :param str out: the path to the merged file
    :param str ignore: a string which start the lines to ignore
    :param str header: The header of the merged file
    """
    def get_header(result: bool, warnings):
        res_or_not = header if result else f"No {header}"
        header_str = f"""# parallel_msf {macsypy.__version__}
# merged {os.path.basename(files[0])}
# {res_or_not} found:
"""
        if warnings:
            header_str += '# \n'
            for warn in warnings:
                header_str += warn
            header_str += '# \n'

        return header_str

    warnings = []
    with open(out, 'w') as f_out:
        last_sol_id = 0
        new_sol_id = None
        first = True
        results = False
        for file in files:
            _log.debug(f"Merging {file}")
            warnings += get_warning(file)
            skip = bool(skip_until)
            with open(file) as fh_in:
                for line in fh_in:
                    if skip:
                        if skip_until(line):
                            skip = False
                        else:
                            continue
                    if first and line.startswith('sol_id'):
                        first = False
                        f_out.write(get_header(True, warnings))
                        new_line = line
                    elif line.startswith('sol_id'):
                        continue
                    elif comment and line.startswith(comment):
                        f_out.write(line)
                        continue
                    elif line.startswith('\n'):
                        new_line = line
                    elif line:
                        fields = line.split('\t')
                        try:
                            new_sol_id = int(fields[0]) + last_sol_id
                        except ValueError as err:
                            msg = f"Cannot reindex int({fields[0]}) + {last_sol_id}: {err}"
                            _log.critical(msg)
                            raise ValueError(msg) from None
                        fields[0] = str(new_sol_id)
                        new_line = '\t'.join(fields)

                    f_out.write(new_line)
                    results = True
            if new_sol_id is not None:  # The previous results are empty
                last_sol_id = new_sol_id
        if not results:
            f_out.write(get_header(False, warnings) + '\n')


def merge_summary(files: List[str], out: str, header: str = "") -> None:
    """

    :param files: the list of files to merge
    :param str out: the path to the merged file
    :param str header: The header of the merged file
    :return:
    """
    warnings = []
    data = []
    for one_file in files:
        datum = pd.read_csv(one_file, sep="\t", comment="#", index_col="replicon")
        data.append(datum)
        warnings += get_warning(one_file)
    merged = pd.concat(data, axis=0)
    with open(out, 'w') as f_out:
        res_or_not = f"# {header}:" if not merged.empty else "# No Systems found:"
        header_str = f"""# parallel_msf {macsypy.__version__}
# merged {os.path.basename(files[0])}
{res_or_not}
"""
        f_out.write(header_str)
        for warn in warnings:
            f_out.write('# \n')
            f_out.write(warn)
            f_out.write('# \n')
        merged.to_csv(f_out, sep="\t")


def merge_results(results_dirs: List[str], out_dir: str = '.') -> None:
    """

    :param results_dirs: The list of macsyfinder results directories to merge
    :type results_dirs: list of str
    :param str out_dir: the path to the directory where to store the merged files
    """
    filename_to_merge, ext = 'all_best_solutions', 'tsv'
    out_file = os.path.join(out_dir, f'merged_{filename_to_merge}.{ext}')
    _log.info(f"Merging '{filename_to_merge}.{ext}' in to '{out_file}'")
    all_best_solutions_files = [os.path.join(d, f'{filename_to_merge}.{ext}') for d in results_dirs]
    header = "Systems"
    merge_and_reindex(all_best_solutions_files, out_file, header,
                      skip_until=lambda l: l.startswith('sol_id'),
                      comment='#')

    for filename_to_merge, ext, header, first_col in [('best_solution', 'tsv', 'Systems', 'replicon'),
                                           ('all_systems', 'tsv', 'Systems', 'replicon'),
                                           ('best_solution_multisystems', 'tsv', 'Multisystems', 'replicon'),
                                           ('best_solution_loners', 'tsv', 'Loners', 'replicon'),
                                           ('rejected_candidates', 'tsv', 'Rejected', 'candidate_id'),
                                           ]:
        out_file = os.path.join(out_dir, f'merged_{filename_to_merge}.{ext}')
        _log.info(f"Merging '{filename_to_merge}.{ext}' in to '{out_file}'")
        best_solution_files = [os.path.join(d, f'{filename_to_merge}.{ext}') for d in results_dirs]
        merge_files(best_solution_files, out_file, header,
                    skip_until=lambda l: l.startswith(first_col),
                    keep_first=first_col)

    filename_to_merge = 'all_systems'
    ext = 'txt'
    out_file = os.path.join(out_dir, f'merged_{filename_to_merge}.{ext}')
    _log.info(f"Merging '{filename_to_merge}.{ext}' in to '{out_file}'")
    filename_to_merge = [os.path.join(d, f'{filename_to_merge}.{ext}') for d in results_dirs]
    merge_files(filename_to_merge, out_file, 'Systems',
                skip_until=lambda l: l.startswith('system id'),
                keep_first='system id')

    filename_to_merge = 'rejected_candidates'
    ext = 'txt'
    out_file = os.path.join(out_dir, f'merged_{filename_to_merge}.{ext}')
    _log.info(f"Merging '{filename_to_merge}.{ext}' in to '{out_file}'")
    filename_to_merge = [os.path.join(d, f'{filename_to_merge}.{ext}') for d in results_dirs]
    merge_files(filename_to_merge, out_file, 'Rejected candidates',
                skip_until=lambda l: l.startswith('Cluster:'),
                keep_first='Cluster:')

    filename_to_merge, ext = 'best_solution_summary', 'tsv'
    out_file = os.path.join(out_dir, f'merged_{filename_to_merge}.{ext}')
    _log.info(f"Merging '{filename_to_merge}' in to '{out_file}'")
    all_summary_files = [os.path.join(d, f'{filename_to_merge}.{ext}') for d in results_dirs]
    merge_summary(all_summary_files, out_file, header='Best Solution Summary')


def parse_args(args:  List[str]) -> argparse.Namespace:
    """

    :param args: the arguments passed on the command line without the first elemnet
    :type args: list of str
    :return: the command line and options arguments parsed
    :rtype: :class:`argparse.Namespace` object
    """
    description = """Merge the different files from several macsyfinder results in one.
    
    - merge the 'best_solution.tsv' in to 'merged_best_solution.tsv'
    - merge the 'best_multisystems.tsv' in to 'merged_best_multisystems.tsv'
    - merge the 'best_loners.tsv' in to 'merged_best_loners.tsv'
    - merge the 'all_best_solutions.tsv' in to `merged_all_best_solutions'
    - merge the 'all_systems.tsv' in to 'merged_all_systems.tsv'
    - merge the 'all_systems.txt' in to 'merged_all_systems.txt'
    - merge the 'rejected_candidates.txt' in to 'merged_rejected_candidates.txt'
"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=description)
    parser.add_argument('results_dirs',
                         nargs='+',
                         help='Path to the macsyfinder results directories to merge eg : path/to/macsyfinder-date-hour')
    parser.add_argument('-o', '--out-dir',
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
    if not os.path.exists(parsed_args.out_dir):
        try:
            os.mkdir(parsed_args.out_dir)
        except PermissionError as err:
            outdir_err = f"Cannot create {parsed_args.out_dir} : {err}"
    elif not os.path.isdir(parsed_args.out_dir):
        outdir_err = f"{parsed_args.out_dir} is not a directory"
    elif not os.access(parsed_args.out_dir, os.W_OK):
        outdir_err = f"{parsed_args.out_dir} is not writable"
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

    macsypy.init_logger(log_file=os.path.join(parsed_args.out_dir, 'macsy_merge_results.out'),
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

    merge_results(parsed_args.results_dirs, out_dir=parsed_args.out_dir)


if __name__ == '__main__':
    main()
