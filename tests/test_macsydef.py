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
import tempfile
import shutil
from xml.etree import ElementTree as Et
from tests import MacsyTest
from macsypy.scripts.macsydef_1to2 import verbosity_to_log_level, parse_args, _1to2, main


class TestMacsydef(MacsyTest):

    def test_verbosity_to_log_level(self):
        self.assertEqual(verbosity_to_log_level(0), 20)
        self.assertEqual(verbosity_to_log_level(1), 10)
        self.assertEqual(verbosity_to_log_level(2), 0)


    def test_parse_args(self):
        args = ['foo.xml', 'bar.xml']
        parsed_args = parse_args(args)
        self.assertListEqual(parsed_args.definitions, ['foo.xml', 'bar.xml'])
        self.assertEqual(parsed_args.verbosity, 0)
        self.assertFalse(parsed_args.in_place)

        args = ['-vv', '-i', 'foo.xml', 'bar.xml']
        parsed_args = parse_args(args)
        self.assertEqual(parsed_args.verbosity, 2)
        self.assertTrue(parsed_args.in_place)


    def test_1to2(self):
        src_xml = self.find_data('models', 'old', 'definitions', 'model_3.xml')
        new_tree = _1to2(src_xml)
        root = new_tree.getroot()
        self.assertEqual(root.tag, 'model')
        sys_ref = root.findall(".//gene[@system_ref]")
        self.assertListEqual(sys_ref, [])
        model_ref = root.findall(".//gene[@model_ref]")
        self.assertEqual(len(model_ref), 3)


    def test_main(self):
        tmpdir = os.path.join(tempfile.gettempdir(), 'test_macsydef')
        os.mkdir(tmpdir)
        src_xml = self.find_data('models', 'old', 'definitions', 'model_4.xml')
        test_xml = shutil.copy(src_xml, tmpdir)
        try:
            main([test_xml], loglevel=30)
            self.assertTrue(os.path.exists(f"{test_xml}.ori"))
            tree = Et.parse(test_xml)
            root = tree.getroot()
            self.assertEqual(root.tag, 'model')
            sys_ref = root.findall(".//gene[@system_ref]")
            self.assertListEqual(sys_ref, [])
            model_ref = root.findall(".//gene[@model_ref]")
            self.assertEqual(len(model_ref), 1)
        finally:
            try:
                shutil.rmtree(tmpdir)
            except Exception:
                pass


    def test_main_in_place(self):
        tmpdir = os.path.join(tempfile.gettempdir(), 'test_macsydef')
        os.mkdir(tmpdir)
        src_xml = self.find_data('models', 'old', 'definitions', 'model_4.xml')
        test_xml = shutil.copy(src_xml, tmpdir)
        try:
            main(['-i', test_xml], loglevel=30)
            self.assertFalse(os.path.exists(f"{test_xml}.ori"))
            tree = Et.parse(test_xml)
            root = tree.getroot()
            self.assertEqual(root.tag, 'model')
            sys_ref = root.findall(".//gene[@system_ref]")
            self.assertListEqual(sys_ref, [])
            model_ref = root.findall(".//gene[@model_ref]")
            self.assertEqual(len(model_ref), 1)
        finally:
            try:
                shutil.rmtree(tmpdir)
            except Exception:
                pass


    def main_with_error(self):
        src_xml = self.find_data(os.path.join('models', 'foo', 'definitions', 'not_xml.xml'))
        with self.catch_log(log_name='macsydef') as log:
            main([src_xml])
            log_msg = log.get_value()
        self.assertEqual(log_msg, f"The definition file {src_xml} "
                                  f"cannot be migrate: syntax error: line 1, column 0 : skip it.")
