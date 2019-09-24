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
import argparse
import xml.etree.ElementTree as Et
import logging
import colorlog
_log = colorlog.getLogger('macsydef')


def init_logger(level='INFO', out=True):
    logger = colorlog.getLogger('macsydef')
    logging = colorlog.logging.logging
    handlers = []
    if out:
        stdout_handler = colorlog.StreamHandler(sys.stderr)
        stdout_formatter = colorlog.ColoredFormatter("%(log_color)s%(message)s",
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
        handlers.append(stdout_handler)
    else:
        null_handler = logging.NullHandler()
        logger.addHandler(null_handler)
        handlers.append(null_handler)
    if isinstance(level, str):
        level = getattr(logging, level)
    logger.setLevel(level)
    return logger


def verbosity_to_log_level(verbosity):
    """
    transform the number of -v option in loglevel
    :param int verbosity: number of -v option on the command line
    :return: an int corresponding to a logging level
    """
    level = max((logging.INFO - (10 * verbosity), 0))
    return level


def parse_args(args):
    """

    :param args: The arguments provided on the command line
    :type args: List of strings [without the program name]
    :return: The arguments parsed
    :rtype: :class:`aprgparse.Nampsace` object.
    """
    parser = argparse.ArgumentParser(description="""Migrate model definition from macsyfinder 1.x syntax to 2.x
By default the old xml file is rename as xx.xml.ori.""",
                                     epilog="For more details, visit the MacSyFinder website and see "
                                            "the MacSyFinder documentation.")
    parser.add_argument('definitions',
                        nargs='+',
                        help="the definitions path to migrate"
                        )
    parser.add_argument('--in-place', '-i',
                        action='store_true',
                        default=False,
                        help="modify xml definition in place.")

    parser.add_argument('-v',
                        dest='verbosity',
                        action="count",
                        default=0,
                        help="Increases the verbosity level."
                        )
    parsed_args = parser.parse_args(args)
    return parsed_args


def _1to2(xml):
    """
    Parse the xml and update the syntax

    :param str xml: the path to the xml file
    :return: the updated tree
    :rtype: :class:`xml.etree.ElementTree` object
    """
    tree = Et.parse(xml)
    root = tree.getroot()
    root.tag = 'model'
    sys_ref = root.findall(".//gene[@system_ref]")
    for gene_node in sys_ref:
        gene_node.attrib['model_ref'] = gene_node.attrib['system_ref']
        del gene_node.attrib['system_ref']
    return tree


def main(args=None, loglevel=None):
    """
    loop over xml definitions in args to update the syntax

    :param args: the arguments passed on the command line without the program name
    :type args: List of string
    :param loglevel: the output verbosity
    :type loglevel: a positive int or a string among 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
    """
    args = sys.argv[1:] if args is None else args
    parsed_args = parse_args(args)

    if loglevel is None:
        loglevel = verbosity_to_log_level(parsed_args.verbosity)
    init_logger(level=loglevel)

    for xml in parsed_args.definitions:
        _log.info(f"migrate {xml}")
        try:
            tree = _1to2(xml)
        except Et.ParseError as err:
            _log.error(f"The definition file {xml} cannot be migrate: {err} : skip it.")
            continue
        if not parsed_args.in_place:
            os.rename(xml, f"{xml}.ori")
        tree.write(xml)
    _log.info("done")


if __name__ == "__main__":
    main()
