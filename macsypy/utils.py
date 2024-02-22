#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2024  Institut Pasteur (Paris) and CNRS.           #
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
Some macsyfinder helper functions
"""

import os
import os.path
import gzip
import contextlib
import subprocess
from itertools import groupby

from .registries import DefinitionLocation, ModelRegistry
from .error import MacsypyError


def get_def_to_detect(models: list[tuple[str, tuple[str]]],
                      model_registry: ModelRegistry) -> tuple[list[DefinitionLocation], str, str]:
    """
    :param models: the list of models to detect as returned by config.models.
    :type models: list of tuple with the following structure:
                  [('model_fqn', ('def1, def2, ...)), ('model_2', ('def1', ...)), ...]
    :param model_registry: the models registry for this run.
    :return: the definitions to parse
    :raise ValueError: if a model name provided in models is not in model_registry.
    """
    root, def_names = models
    root = root.rstrip(os.path.sep)
    model_family = DefinitionLocation.root_name(root)
    model_loc = model_registry[model_family]
    model_vers = model_loc.version
    if 'all' in [d.lower() for d in def_names]:
        if root == model_loc.name:
            root = None
        def_to_detect = model_loc.get_all_definitions(root_def_name=root)
    else:
        def_to_detect = [model_loc.get_definition(f'{root}/{one_def}') for one_def in def_names]
    return def_to_detect, model_family, model_vers


def get_replicon_names(genome_path, db_type) -> list[str]:
    if db_type == 'gembase':
        return _get_gembase_replicon_names(genome_path)
    elif db_type in ('ordered_replicon', 'unordered'):
        return [os.path.splitext(os.path.basename(genome_path))[0]]
    else:
        raise MacsypyError(f"Invalid genome type: {db_type}")


def _get_gembase_replicon_names(genome_path: str) -> list[str]:
    """
    parse gembase file and get the list of replicon identifiers

    :param genome_path: The path to a file containing sequence in **gembase** format
    :return: the list of replicon identifiers
    """
    def grp_replicon(ids: str) -> str:
        """
        in gembase the identifier of fasta sequence follows the following schema:
        <replicon-name>_<seq-name> with eventually '_' inside the <replicon_name>
        but not in the <seq-name>.
        so grp_replicon allow to group sequences belonging to the same replicon.
        """
        return "_".join(ids.split('_')[: -1])

    seq_ids = []
    with open(genome_path, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                seq_ids.append(line.split()[0][1:])
    replicons = [rep_name for rep_name, _ in groupby(seq_ids, key=grp_replicon)]
    return replicons


def threads_available() -> int:
    """

    :return: The maximal number of threads available.
             It's nice with cluster scheduler or linux.
             On Mac it uses the number of physical cores
    """
    if hasattr(os, "sched_getaffinity"):
        threads_nb = len(os.sched_getaffinity(0))
    else:
        threads_nb = os.cpu_count()
    return threads_nb


def parse_time(user_time: int | str) -> int:
    """
    parse user-friendly time and return it in seconds
    user time supports units as s h m d for sec min hour day
    or a combination of them
    1h10m50s means 1 hour 10 minutes 50 seconds
    all terms will be converted in seconds and added

    :param user_time:
    :return: seconds
    :raise: ValueError if user_time is not parseable
    """
    try:
        user_time = int(user_time)
        return user_time  # user time has no units , it's seconds
    except ValueError:
        pass
    import re
    parts_converter = {'s': lambda x: x,
                       'm': lambda x: x * 60,
                       'h': lambda x: x * 3600,
                       'd': lambda x: x * 86400
                       }
    time_parts = re.findall(r'(\d+)(\D+)', user_time)
    time = 0
    for value, unit in time_parts:
        unit = unit.strip().lower()
        try:
            time += parts_converter[unit](int(value))
        except KeyError:
            raise ValueError("Not valid time format. Units allowed h/m/s.")
    return time


@contextlib.contextmanager
def open_compressed(path: str, mode: str = 'rt') -> str:
    """

    :param path: the path to open
    :param mode: the opening mode by default read text
    :yield: the content of the file line by line
    """
    _, ext = os.path.splitext(path)
    if ext == '.gz':
        my_open = gzip.open
    elif ext == '.bz2' or ext == '.zip':
        msg = f"MacSyFinder does not support '{ext[1:]}' compression (only gzip)."
        raise ValueError(msg)
    else:
        # I assumed it's a fasta not compressed
        my_open = open
    with my_open(path, mode) as f:
        yield f


def get_git_revision_short_hash():
    """
    :return: the git commit number (short version)
    :rtype: str
    """
    try:
        short_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'],
                                             cwd=os.path.dirname(os.path.abspath(__file__)))
    except subprocess.CalledProcessError:
        short_hash = 'not_git'
    short_hash = str(short_hash, "utf-8").strip()
    return short_hash
