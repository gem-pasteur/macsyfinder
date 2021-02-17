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

import contextlib
import os
import sys

from .registries import DefinitionLocation


def get_def_to_detect(models, model_registry):
    """
    :param models: the list of models to detect as returned by config.models.
    :type models: list of tuple with the following structure:
                  [('model_1', ('def1, def2, ...)), ('model_2', ('def1', ...)), ...]
    :param model_registry: the models registry for this run.
    :type model_registry: :class:`macsypy.registries.ModelRegistry` object.
    :return: the definitions to parse
    :rtype: list of :class:`macsypy.registries.DefinitionLocation` objects
    :raise ValueError: if a model name provided in models is not in model_registry.
    """
    root, def_names = models
    model_family = DefinitionLocation.root_name(root)
    model_loc = model_registry[model_family]
    if 'all' in [d.lower() for d in def_names]:
        if root == model_loc.name:
            root = None
        def_to_detect = model_loc.get_all_definitions(root_def_name=root)
    else:
        def_to_detect = [model_loc.get_definition(f'{root}/{one_def}') for one_def in def_names]
    return def_to_detect

