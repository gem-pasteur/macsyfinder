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

import xml.etree.ElementTree as Et
import logging
_log = logging.getLogger(__name__)

from macsypy.error import MacsypyError


class ModelConfParser:
    """
    Handle model_conf.xml configuration file.
    """

    def __init__(self, path):
        """

        :param str path: The path to the configuration file
        """
        self._path = path


    def parse(self):
        """
        Parse the xml 'model_conf' file set at the root of a data package

        :return: The specific configuration for a model family
        :rtype: dict with the name of variables as keys and value as values
        """
        model_conf_node = self._get_model_conf_node()
        weights_node = model_conf_node.find("./weights")

        filtering_opt = {}
        weights = {}
        if weights_node:
            weights = self.parse_weights(weights_node)

        filtering_node = model_conf_node.find("./filtering")
        if filtering_node:
            filtering_opt = self.parse_filtering(filtering_node)

        model_conf = {k: v for conf_part in (weights, filtering_opt) for k, v in conf_part.items()}
        return model_conf


    def _get_model_conf_node(self):
        """
        Find the root of the document

        :return: the document root of model_conf
        """
        try:
            tree = Et.parse(self._path)
            model_node = tree.getroot()
        except Exception as err:
            msg = f"unable to parse model configuration '{self._path}' : {err}"
            _log.critical(msg)
            raise MacsypyError(msg) from None
        return model_node


    def parse_weights(self, weights_node):
        """
        Parse the node 'weights' contening the scoring weight configuration

        :param weights_node: the node 'weights'
        :type weights_node: :class"`Et.ElementTree` object
        :return: the configuration option/value about the scores
        :rtype: dict
        """
        elements = {'itself': float,
                    'exchangeable': float,
                    'mandatory': float,
                    'accessory': float,
                    'neutral': float,
                    'out_of_cluster': float,
                    'redundancy_penalty': float}
        weights_conf = self._parse_section(weights_node, elements)
        # rename options as in the other part of MSF
        weights_conf = {(f"{k}_weight"if k != 'redundancy_penalty' else k): v for k, v in weights_conf.items()}
        return weights_conf


    def parse_filtering(self, filtering_node):
        """
        Parse the node 'filtering' containing the filtering options configuration

        :param filtering_node: the node 'filtering'
        :type filtering_node: :class"`Et.ElementTree` object
        :return: the configuration option/value about the filtering
        :rtype: dict
        """
        def parse_cut_ga(value):
            if value.lower() in ('true', 1):
                return True
            elif value.lower() in ('false', 0):
                return False
            else:
                msg = f"cannot parse 'cut_ga' element in '{self._path}' expect True, 1, False, 0 got : '{value}'"
                _log.critical(msg)
                raise MacsypyError(msg)

        elements = {'e_value_search': float,
                    'i_evalue_sel': float,
                    'coverage_profile': float,
                    'cut_ga': parse_cut_ga,
                    }
        fiter_conf = self._parse_section(filtering_node, elements)
        return fiter_conf


    def _parse_section(self, section_node, allowed_elements):
        """
        Parse a node containing configurations options and value

        :param section_node:
        :param allowed_elements: The elements allowed in this section
                                 Only these elements are parsed and in the final dictionnary
        :type allowed_elements: a dict with options name as keys and function to parse the element
        :return: dict
        """
        section = {}
        for child in section_node:
            element = child.tag
            if element in allowed_elements:
                value = child.text
                try:
                    value = allowed_elements[element](value)
                except (TypeError, ValueError) as err:
                    msg = f"The model configuration file '{self._path}' cannot be parsed: {err}"
                    _log.critical(msg)
                    raise MacsypyError(msg) from None
            else:
                _log.warning(f"unknown element '{element}' in '{self._path}' ignore it.")
                continue
            section[element] = value
        return section
