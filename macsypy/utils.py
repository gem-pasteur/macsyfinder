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
Some macsyfinder helper functions
"""

import os
import os.path
from itertools import groupby

from .registries import DefinitionLocation
from .error import MacsypyError


def get_def_to_detect(models, model_registry):
    """
    :param models: the list of models to detect as returned by config.models.
    :type models: list of tuple with the following structure:
                  [('model_fqn', ('def1, def2, ...)), ('model_2', ('def1', ...)), ...]
    :param model_registry: the models registry for this run.
    :type model_registry: :class:`macsypy.registries.ModelRegistry` object.
    :return: the definitions to parse
    :rtype: list of :class:`macsypy.registries.DefinitionLocation` objects
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


def get_replicon_names(genomee_path, db_type):
    if db_type == 'gembase':
        return _get_gembase_replicon_names(genomee_path)
    elif db_type in ('ordered_replicon', 'unordered'):
        return [os.path.splitext(os.path.basename(genomee_path))[0]]
    else:
        raise MacsypyError(f"Invalid genome type: {db_type}")


def _get_gembase_replicon_names(genome_path):
    """
    parse gembase file and get the list of replicon identifiers

    :param str genome_path: The path to a file containing sequence in **gembase** format
    :return: the list of replicon identifiers
    :rtype: list of str
    """
    def grp_replicon(ids):
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


def threads_available():
    """

    :return: The maximal number of threads available.
             It's nice with cluster scheduler or linux.
             On Mac it use the number of physical cores
    :rtype: int
    """
    if hasattr(os, "sched_getaffinity"):
        threads_nb = len(os.sched_getaffinity(0))
    else:
        threads_nb = os.cpu_count()
    return threads_nb


def indent_wrapper(ElementTree):
    """
    xml.etree.ElementTree implement ident only from python 3.9
    below the code from python 3.9 to inject it in ET at runtime

    :param ElementTree: ElementTree class
    :type ElementTree: class
    :return: function indent
    :rtype: function
    """

    def indent(tree, space="  ", level=0):
        """Indent an XML document by inserting newlines and indentation space
        after elements.

        *tree* is the ElementTree or Element to modify.  The (root) element
        itself will not be changed, but the tail text of all elements in its
        subtree will be adapted.

        *space* is the whitespace to insert for each indentation level, two
        space characters by default.

        *level* is the initial indentation level. Setting this to a higher
        value than 0 can be used for indenting subtrees that are more deeply
        nested inside of a document.
        """
        if isinstance(tree, ElementTree):
            tree = tree.getroot()
        if level < 0:
            raise ValueError(f"Initial indentation level must be >= 0, got {level}")
        if not len(tree):
            return

        # Reduce the memory consumption by reusing indentation strings.
        indentations = ["\n" + level * space]

        def _indent_children(elem, level):
            # Start a new indentation level for the first child.
            child_level = level + 1
            try:
                child_indentation = indentations[child_level]
            except IndexError:
                child_indentation = indentations[level] + space
                indentations.append(child_indentation)

            if not elem.text or not elem.text.strip():
                elem.text = child_indentation

            for child in elem:
                if len(child):
                    _indent_children(child, child_level)
                if not child.tail or not child.tail.strip():
                    child.tail = child_indentation

            # Dedent after the last child by overwriting the previous indentation.
            if not child.tail.strip():
                child.tail = indentations[level]

        _indent_children(tree, 0)

    return indent


def parse_time(user_time):
    """
    parse user friendly time and return it in seconds
    user time supports units as s h m d for sec min hour day
    or a combination of them
    1h10m50s means 1 hour 10 minutes 50 seconds
    all terms will be converted in seconds and added

    :param user_time:
    :type user_time: int or str
    :return: seconds
    :rtype: int
    :raise: ValueError if user_time is not parseable
    """
    try:
        user_time = int(user_time)
        return user_time # user time has no units , it's seconds
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

