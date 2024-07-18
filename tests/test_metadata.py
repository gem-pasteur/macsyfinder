##########################################################################
#  MacSyFinder - Detection of macromolecular systems in protein dataset  #
#                using systems modelling and similarity search.          #
#  Authors: Sophie Abby, Bertrand Neron                                  #
#  Copyright (c) 2014-2024  Institut Pasteur (Paris) and CNRS.           #
#  See the COPYRIGHT file for details                                    #
#                                                                        #
#  This file is part of MacSyFinder package.                             #
#                                                                        #
#  MacSyFinder is free software: you can redistribute it and/or modify   #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  MacSyFinder is distributed in the hope that it will be useful,        #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details .                         #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with MacSyFinder (COPYING).                                     #
#  If not, see <https://www.gnu.org/licenses/>.                          #
##########################################################################

import os
import shutil
import tempfile

import yaml

from macsypy.metadata import Metadata, Maintainer
from tests import MacsyTest


class TestMetadata(MacsyTest):

    def setUp(self) -> None:
        self.tmpdir = os.path.join(tempfile.gettempdir(), 'macsy_test_metadata')
        if os.path.exists(self.tmpdir) and os.path.isdir(self.tmpdir):
            shutil.rmtree(self.tmpdir)
        os.makedirs(self.tmpdir)


    def tearDown(self) -> None:
        try:
            shutil.rmtree(self.tmpdir)
        except Exception:
            pass


    def test_init(self):
        maintainer = Maintainer(name='joe', email='joe@bar.org')
        short_desc = 'this is a short description'
        meta = Metadata(maintainer, short_desc)

        self.assertEqual(meta.maintainer, maintainer)
        self.assertEqual(meta.short_desc, short_desc)
        self.assertIsNone(meta.vers)
        self.assertEqual(meta.license, '')


    def test_load(self):
        meta_path = self.find_data('pack_metadata', 'good_metadata.yml')
        meta = Metadata.load(meta_path)
        self.assertEqual(meta.maintainer, Maintainer(name='auth_name', email='auth_name@mondomain.fr'))
        self.assertEqual(meta.short_desc, 'this is a short description of the repos')
        self.assertEqual(meta.vers, '0.0b2')
        self.assertEqual(meta.license, 'CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)')
        self.assertEqual(meta.doc, 'http://link/to/the/documentation')
        self.assertEqual(meta.cite, ['bla bla',
                                                'link to publication',
                                                'ligne 1\nligne 2\nligne 3 et bbbbb\n'])

        meta_path = self.find_data('pack_metadata', 'bad_metadata.yml')
        with self.assertRaises(ValueError) as ctx:
            Metadata.load(meta_path)
        self.assertEqual(str(ctx.exception),
                         f"""
- The metadata file '{meta_path}' is not valid: the element 'short_desc' is required.
- The metadata file '{meta_path}' is not valid: the element 'maintainer' is required.""")


    def test_maintainer(self):
        meta_path = self.find_data('pack_metadata', 'good_metadata.yml')
        meta = Metadata.load(meta_path)
        self.assertEqual(meta.maintainer, Maintainer(name='auth_name', email='auth_name@mondomain.fr'))
        new_auth = 'joe bar'
        new_email = 'joe@domain.org'
        meta.maintainer = Maintainer(new_auth, new_email)
        self.assertEqual(meta.maintainer, Maintainer(name=new_auth, email=new_email))


    def test_short_desc(self):
        meta_path = self.find_data('pack_metadata', 'good_metadata.yml')
        meta = Metadata.load(meta_path)
        self.assertEqual(meta.short_desc, 'this is a short description of the repos')
        new_short_desc = "This is a new short desc"
        meta.short_desc = new_short_desc
        self.assertEqual(meta.short_desc, new_short_desc)
        with self.assertRaises(ValueError) as ctx:
            meta.short_desc = ''
        self.assertEqual(str(ctx.exception),
                         "The field 'short_desc' is mandatory.")


    def test_vers(self):
        meta_path = self.find_data('pack_metadata', 'good_metadata.yml')
        meta = Metadata.load(meta_path)
        self.assertEqual(meta.vers, '0.0b2')
        meta.vers = '0.0b3'
        self.assertEqual(meta.vers, '0.0b3')


    def test_cite(self):
        meta_path = self.find_data('pack_metadata', 'good_metadata.yml')
        meta = Metadata.load(meta_path)
        self.assertEqual(meta.cite, ['bla bla', 'link to publication', 'ligne 1\nligne 2\nligne 3 et bbbbb\n'])
        new_cite = ["citation_1", "citation_2"]
        meta.cite = new_cite
        self.assertEqual(meta.cite, new_cite)


    def test_license(self):
        meta_path = self.find_data('pack_metadata', 'good_metadata.yml')
        meta = Metadata.load(meta_path)
        self.assertEqual(meta.license, 'CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)')
        new_license = "GPLv3"
        meta.license = new_license
        self.assertEqual(meta.license, new_license)


    def test_copyright(self):
        meta_path = self.find_data('pack_metadata', 'good_metadata.yml')
        meta = Metadata.load(meta_path)
        self.assertEqual(meta.copyright, '2019, Institut Pasteur, CNRS')
        self.assertEqual(meta.copyright_date, '2019')
        self.assertEqual(meta.copyright_holder, 'Institut Pasteur, CNRS')
        meta.copyright_date = '2020-2024'
        meta.copyright_holder = 'My institution'
        self.assertEqual(meta.copyright, '2020-2024, My institution')


    def test_doc(self):
        meta_path = self.find_data('pack_metadata', 'good_metadata.yml')
        meta = Metadata.load(meta_path)
        self.assertEqual(meta.doc, 'http://link/to/the/documentation')
        meta.doc = 'new link'
        self.assertEqual(meta.doc, 'new link')


    def test_save(self):
        for yml_name in ('good_metadata.yml', 'metadata_no_vers.yml', 'metadata_no_license.yml',
                         'metadata_no_doc.yml', 'metadata_no_copyright.yml'):
            with self.subTest(yml_name=yml_name):
                meta_path = self.find_data('pack_metadata', yml_name)
                meta = Metadata.load(meta_path)
                saved_path = os.path.join(self.tmpdir, 'metadata.yml')
                meta.save(saved_path)
                with open(saved_path) as saved_file:
                    saved_yaml = yaml.safe_load(saved_file)
                with open(meta_path) as ref_file:
                    ref_yaml = yaml.safe_load(ref_file)
                self.assertDictEqual(ref_yaml, saved_yaml)
