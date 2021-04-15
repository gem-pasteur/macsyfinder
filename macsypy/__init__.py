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

from time import strftime, localtime
import sys

__version__ = "2.0rc2"

__MACSY_CONF__ = '$MACSYCONF'
__MACSY_DATA__ = '$MACSYDATA'

__citation__ = """Abby SS, Néron B, Ménager H, Touchon M, Rocha EPC (2014)
MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems.
PLoS ONE 9(10): e110726. doi:10.1371/journal.pone.0110726
"""


def init_logger(log_file=None, out=True):
    """

    :param str log_file: The path toward a file log
    :param out:
    :return:
    """
    import colorlog
    logger = colorlog.getLogger('macsypy')
    logging = colorlog.logging.logging
    handlers = []
    if out:
        stdout_handler = colorlog.StreamHandler(sys.stdout)
        stdout_formatter = colorlog.ColoredFormatter("%(log_color)s%(message)s",
                                                     datefmt=None,
                                                     reset=True,
                                                     log_colors={
                                                         'DEBUG':    'cyan',
                                                         'INFO':     'green',
                                                         'WARNING':  'yellow',
                                                         'ERROR':    'red',
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
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_formatter = logging.Formatter("%(levelname)-8s : %(message)s")
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        handlers.append(file_handler)
    logger.setLevel(logging.WARNING)
    return handlers


def logger_set_level(level='INFO'):
    # default value must be a string
    # cannot be colorlog.logging.logging.WARNING for instance
    # because setup import __init__ to get __version__
    # so logger_set_level is defined
    # if level is colorlog.logging.logging.WARNING
    # that mean that colorlog must be already installed
    # otherwise an error occured during pip install
    #  NameError: name 'colorlog' is not defined
    import colorlog

    levels = {'NOTSET': colorlog.logging.logging.NOTSET,
              'DEBUG': colorlog.logging.logging.DEBUG,
              'INFO': colorlog.logging.logging.INFO,
              'WARNING': colorlog.logging.logging.WARNING,
              'ERROR': colorlog.logging.logging.ERROR,
              'CRITICAL': colorlog.logging.logging.CRITICAL,
              }
    if level in levels:
        level = levels[level]
    elif not isinstance(level, int) and level < 0:
        raise ValueError(f"Level must be {', '.join(levels.keys())} or a positive integer")

    logger = colorlog.getLogger('macsypy')
    if level <= colorlog.logging.logging.DEBUG:
        stdout_formatter = colorlog.ColoredFormatter(
            "%(log_color)s%(levelname)-8s : %(module)s: L %(lineno)d :%(reset)s %(message)s",
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
        stdout_handler = logger.handlers[0]
        stdout_handler.setFormatter(stdout_formatter)

        logging = colorlog.logging.logging
        if len(logger.handlers) > 1:
            file_formatter = logging.Formatter("%(levelname)-8s : %(module)s: L %(lineno)d : %(message)s")
            file_handler = logger.handlers[1]
            file_handler.setFormatter(file_formatter)

    logger.setLevel(level)
