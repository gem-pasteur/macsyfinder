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

import logging
import colorlog
import xml.etree.ElementTree as Et

import macsypy
from macsypy.error import MacsypyError
from macsypy.model_conf_parser import ModelConfParser
import macsypy.model_conf_parser
from tests import MacsyTest


class TestModelConfParser(MacsyTest):

    def setUp(self) -> None:
        # need to do hugly trick with logger
        # because logger are singleton and
        # trigger some side effect with othe unit tests
        # for instance if run the test below after test_macsypy
        # where I tests loggers model_conf_parser
        # is not replaced by the logger set in setup
        # then the catch_log doesn't work anymore
        macsypy.init_logger()
        macsypy.logger_set_level(logging.INFO)
        logger = colorlog.getLogger('macsypy')
        macsypy.model_conf_parser._log = logger

    def tearDown(self) -> None:
        logger = colorlog.getLogger('macsypy')
        del logger.manager.loggerDict['macsypy']


    def test_parse(self):
        expected_conf = {'itself_weight': 11.0,
                         'exchangeable_weight': 12.0,
                         'mandatory_weight': 13.0,
                         'accessory_weight': 14.0,
                         'neutral_weight': 0.0,
                         'out_of_cluster_weight': 10.0,
                         'redundancy_penalty': 20.0,
                         'e_value_search': 0.12,
                         'i_evalue_sel': 0.012,
                         'coverage_profile': 0.55,
                         'cut_ga': False}
        conf_file = self.find_data('conf_files', 'model_conf.xml')
        mcp = ModelConfParser(conf_file)
        test_conf = mcp.parse()
        self.assertDictEqual(expected_conf, test_conf)


    def test_parse_wo_filtering(self):
        expected_conf = {'itself_weight': 11.0,
                         'exchangeable_weight': 12.0,
                         'mandatory_weight': 13.0,
                         'accessory_weight': 14.0,
                         'neutral_weight': 0.0,
                         'redundancy_penalty': 20.0,
                         'out_of_cluster_weight': 10.0}

        conf_file = self.find_data('conf_files', 'model_conf_wo_filtering.xml')
        mcp = ModelConfParser(conf_file)
        test_conf = mcp.parse()
        self.assertDictEqual(expected_conf, test_conf)


    def test_parse_wo_weights(self):
        expected_conf = {'e_value_search': 0.12,
                         'i_evalue_sel': 0.012,
                         'coverage_profile': 0.55,
                         'cut_ga': False}
        conf_file = self.find_data('conf_files', 'model_conf_wo_weights.xml')
        mcp = ModelConfParser(conf_file)
        test_conf = mcp.parse()
        self.assertDictEqual(expected_conf, test_conf)


    def test_parse_bad_file(self):
        conf_file = self.find_data('conf_files', 'project.conf')
        mcp = ModelConfParser(conf_file)

        with self.catch_log(log_name='macsypy'):
            with self.assertRaises(MacsypyError) as ctx:
                mcp.parse()

        self.assertEqual(str(ctx.exception),
                         f"unable to parse model configuration '{conf_file}' : syntax error: line 2, column 0"
                         )


    def test_parse_weights(self):
        expected_weights = {'itself_weight': 11.0,
                            'exchangeable_weight': 12.0,
                            'mandatory_weight': 13.0,
                            'accessory_weight': 14.0,
                            'neutral_weight': 0.0,
                            'redundancy_penalty': 20.0,
                            'out_of_cluster_weight': 10.0,
                            }
        conf_file = self.find_data('conf_files', 'model_conf.xml')
        tree = Et.parse(conf_file)
        model_node = tree.getroot()
        mcp = ModelConfParser(conf_file)
        recieved_weights = mcp.parse_weights(model_node.find("./weights"))
        self.assertDictEqual(expected_weights, recieved_weights)


    def test_parse_filtering(self):
        expected_filters = {'e_value_search': 0.12,
                            'i_evalue_sel': 0.012,
                            'coverage_profile': 0.55,
                            'cut_ga': False
                            }
        conf_file = self.find_data('conf_files', 'model_conf.xml')
        tree = Et.parse(conf_file)
        model_node = tree.getroot()
        mcp = ModelConfParser(conf_file)
        filtering_node = model_node.find("./filtering")
        recieved_filters = mcp.parse_filtering(filtering_node)
        self.assertDictEqual(recieved_filters, expected_filters)

        # test cut_ga True
        cut_ga_node = filtering_node.find('cut_ga')
        cut_ga_node.text = 'True'
        recieved_filters = mcp.parse_filtering(filtering_node)
        self.assertTrue('cut_ga' in recieved_filters)

        # test bad value for cut_ga
        cut_ga_node.text = 'FOO'
        with self.catch_log(log_name='macsypy') as log:
            with self.assertRaises(MacsypyError) as ctx:
                mcp.parse_filtering(filtering_node)

        self.assertEqual(str(ctx.exception),
                         f"cannot parse 'cut_ga' element in '{conf_file}' expect True, 1, False, 0 got : 'FOO'"
                         )

        # test no cut_ga element
        filtering_node.remove(cut_ga_node)
        recieved_filters = mcp.parse_filtering(filtering_node)
        del(expected_filters['cut_ga'])
        self.assertDictEqual(expected_filters, recieved_filters)


    def test_parse_w_unkown_tag(self):
        expected_filters = {'e_value_search': 0.12,
                            'i_evalue_sel': 0.012,
                            'coverage_profile': 0.55,
                            'cut_ga': False
                            }
        conf_file = self.find_data('conf_files', 'model_conf.xml')
        tree = Et.parse(conf_file)
        model_node = tree.getroot()
        mcp = ModelConfParser(conf_file)
        filter_node = model_node.find("./filtering")
        extra_elt = Et.SubElement(filter_node, 'nimportnaoik')
        extra_elt.text = 'Foo'
        with self.catch_log(log_name='macsypy') as log:
            recieved_filters = mcp.parse_filtering(model_node.find("./filtering"))
            log_msg = log.get_value().strip()

        self.assertEqual(log_msg,
                         f"unknown element 'nimportnaoik' in '{conf_file}' ignore it."
                         )

        self.assertDictEqual(recieved_filters, expected_filters)


    def test_parse_w_bad_value(self):
        conf_file = self.find_data('conf_files', 'model_conf.xml')
        tree = Et.parse(conf_file)
        model_node = tree.getroot()
        mcp = ModelConfParser(conf_file)
        filter_node = model_node.find("./filtering")
        coverage_node = filter_node.find('coverage_profile')
        coverage_node.text = "FOO"

        with self.catch_log(log_name='macsypy') as log:
            with self.assertRaises(MacsypyError) as ctx:
                mcp.parse_filtering(model_node.find("./filtering"))
        self.assertEqual(str(ctx.exception),
                         f"The model configuration file '{conf_file}' cannot be parsed: "
                         f"could not convert string to float: 'FOO'"
                         )