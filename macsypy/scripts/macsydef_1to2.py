################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################

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


def verbosity_to_log_level(verbosity: int) -> str:
    level = max((logging.INFO - (10 * verbosity), 0))
    return level


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('definitions',
                        nargs='+',
                        help="the definitions path to migrate"
                        )
    parser.add_argument('-v',
                        dest='verbosity',
                        action="count",
                        default=0,
                        help="Increases the verbosity level."
                        )
    parsed_args = parser.parse_args(args)
    return parsed_args


def _2to3(xml):
    tree = Et.parse(xml)
    root = tree.getroot()
    root.tag = 'model'
    sys_ref = root.findall(".//gene[@system_ref]")
    for gene_node in sys_ref:
        gene_node.attrib['model_ref'] = gene_node.attrib['system_ref']
        del gene_node.attrib['system_ref']
    return tree


def main(args=None, loglevel=None):
    args = sys.argv[1:] if args is None else args
    parsed_args = parse_args(args)

    if loglevel is None:
        loglevel = verbosity_to_log_level(parsed_args.verbosity)
    init_logger(level=loglevel)

    for xml in parsed_args.definitions:
        _log.info(f"migrate {xml}")
        tree = _2to3(xml)
        os.rename(xml, f"{xml}.ori")
        tree.write(xml)
    _log.info("done")


if __name__ == "__main__":
    main()