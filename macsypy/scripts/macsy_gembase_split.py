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

import sys
import os.path
import argparse
from itertools import groupby

import colorlog
from collections import OrderedDict


import macsypy
from macsypy.config import MacsyDefaults


def copy_chunk(fh_in, out, start, stop):
    """
    Copy file from fh_in to out from position start to stop

    :param fh_in: the source file
    :type fh_in: file like object
    :param str out: the destination file name
    :param int start: the position to start the copy
    :param stop: the position to end the copy
    """
    chunk_size = 1024
    fh_in.seek(start)
    with open(out, 'w') as f_out:
        for start in range(start, stop + 1, chunk_size):
            read_size = min((stop - start), chunk_size)
            content = fh_in.read(read_size)
            f_out.write(content)


def split(seq_index, genome_path, outdir='.'):
    """
    split a file with different replicons in gembase format
    in several files with one replicon per file

    :seq_index: the sequences index
    :type seq_index: dict [str seq_id] : (int start, int stop)
    :param str genome_path: the path to the file to split
    :param str outdir: the path of the directory where to write the replicons
    :return: the list of created replicons files
    :rtype: list of string
    """
    def grp_replicon(line):
        """
        in gembase the identifier of fasta sequence follows the following schema:
        <replicon-name>_<seq-name> with eventually '_' inside the <replicon_name>
        but not in the <seq-name>.
        so grp_replicon allow to group sequences belonging to the same replicon.
        """
        return "_".join(line.split('_')[: -1])


    all_seq_files = []
    with open(genome_path, 'r') as fh_in:
        for rep_name, seq_ids in groupby(seq_index.keys(), key=grp_replicon):
            seqs_ids = [id_ for id_ in seq_ids]
            seq_file = os.path.normpath(os.path.join(outdir, f"{rep_name}.fasta"))
            start = seq_index[seqs_ids[0]][0]
            stop = seq_index[seqs_ids[-1]][1]
            _log.info(f"Writing replicon {seq_file}")
            copy_chunk(fh_in, seq_file, start, stop)
            all_seq_files.append(seq_file)
    return all_seq_files


def index_seq(genome_path):
    """
    Index the sequence in the file represented by genome_path

    :param str genome_path: the path to a file containing several sequences in fasta format
    :return: the sequences index
    :rtype: dict [str seq_id] : (int start, int stop)
    """
    index = OrderedDict()
    with open(genome_path, 'r') as fh:
        start = None
        end = None
        line = fh.readline()
        while line:
            if line.startswith('>') and start is not None:
                end = fh.tell() - len(line)
                index[_id] = start, end
                start = end
                _id = line.split()[0][1:]
            elif line.startswith('>'):  # and start is None
                # The first sequence
                start = fh.tell() - len(line)
                _id = line.split()[0][1:]
            line = fh.readline()
        # this is the end of file
        # add the last sequence
        end = fh.tell()
        index[_id] = start, end

    return index


def parse_args(args):
    """
    :param args: The arguments passed on the command line (without the name of the program)
                 Typically sys.argv[1:]
    :type args: list of string.
    :return: the arguments parsed.
    :rtype: a :class:`argparse.Namespace` object.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="Split a gembase protein file in several files, one per replicon." )
    parser.add_argument('genome_path',
                        help='Path to the genomes file (in gembase format), eg : path/to/file.fst or file.fst')

    parser.add_argument('-o', '--outdir',
                        default='.',
                        help='The path to the directory where to write the chunks.')
    parser.add_argument("--mute",
                        action='store_true',
                        default=False,
                        help="mute the log on stdout."
                             "(continue to log on macsy_gembase_split.out)")

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


def main(args=None, log_level=None):
    """
    main entry point to macsy_gembase_split

        1. index the gembase file to identify start/end of each replicon
        2. use this information to split gembase in several files one per replicon
     
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
        log = colorlog.getLogger('macsypy.split')
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

    macsypy.init_logger(log_file=os.path.join(parsed_args.outdir, 'macsy_gembase_split.out'),
                        out=not parsed_args.mute)
    _log = colorlog.getLogger('macsypy.split')

    if not log_level:
        # logs are specify from args options
        config = MacsyDefaults()
        log_level = max(config.log_level - (10 * parsed_args.verbose) + (10 * parsed_args.quiet), 1)
        macsypy.logger_set_level(log_level)
    else:
        # used by unit tests to mute or unmute logs
        macsypy.logger_set_level(log_level)

    idx = index_seq(parsed_args.genome_path)
    replicon_names = split(idx, parsed_args.genome_path, outdir=parsed_args.outdir)
    print(' '.join(replicon_names))


if __name__ == '__main__':
    main()
