# -*- coding: utf-8 -*-

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

import tempfile
import shutil
import os
import argparse

from macsypy.scripts import macsydata

from tests import MacsyTest


class TestMacsydata(MacsyTest):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.args = argparse.Namespace()
        self.args.org = 'foo'
        self._remote_exists = macsydata.RemoteModelIndex.remote_exists
        macsydata.RemoteModelIndex.remote_exists = lambda x: True

    def tearDown(self):
        macsydata.RemoteModelIndex.remote_exists = self._remote_exists
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass

    def create_fake_package(self, model,
                            xml=True,
                            hmm=True,
                            metadata=True,
                            readme=True,
                            licence=True):
        pack_path = os.path.join(self.tmpdir, model)
        os.mkdir(pack_path)
        if xml:
            def_dir = os.path.join(pack_path, 'definitions')
            os.mkdir(def_dir)
            with open(os.path.join(def_dir, "model_1.xml"), 'w') as f:
                f.write("""<system inter_gene_max_space="20" min_mandatory_genes_required="1" min_genes_required="2">
    <gene name="flgB" presence="mandatory"/>
    <gene name="flgC" presence="mandatory" inter_gene_max_space="2"/>
</system>""")
            with open(os.path.join(def_dir, "model_2.xml"), 'w') as f:
                f.write("""<system inter_gene_max_space="20" min_mandatory_genes_required="1" min_genes_required="2">
    <gene name="fliE" presence="mandatory" multi_system="True"/>
    <gene name="tadZ" presence="accessory" loner="True"/>
    <gene name="sctC" presence="forbidden"/>
</system>""")

        if hmm:
            profile_dir = os.path.join(pack_path, 'profiles')
            os.mkdir(profile_dir)
            for name in ('flgB', 'flgC', 'fliE', 'tadZ', 'sctC'):
                open(os.path.join(profile_dir, f"{name}.hmm"), 'w').close()
        if metadata:
            meta_file = self.find_data('pack_metadata', 'good_metadata.yml')
            meta_dest = os.path.join(pack_path, 'metadata.yml')
            shutil.copyfile(meta_file, meta_dest)
        if readme:
            with open(os.path.join(pack_path, "README"), 'w') as f:
                f.write("# This a README\n")
        if licence:
            with open(os.path.join(pack_path, "LICENCE"), 'w') as f:
                f.write("# This the Licence\n")
        return pack_path


    def test_available(self):
        list_pack = macsydata.RemoteModelIndex.list_packages
        list_pack_vers = macsydata.RemoteModelIndex.list_package_vers
        meta = macsydata.RemoteModelIndex.get_metadata
        pack_name = 'fake_model'
        pack_vers = '1.0'
        pack_meta = {'short_desc': 'desc about fake_model'}
        macsydata.RemoteModelIndex.list_packages = lambda x: [pack_name]
        macsydata.RemoteModelIndex.list_package_vers = lambda x, pack: pack_vers
        macsydata.RemoteModelIndex.get_metadata = lambda x, pack, vers: pack_meta
        self.create_fake_package('fake_model')
        try:
            with self.catch_io(out=True) as out:
                macsydata.do_available(self.args)
                get_pack = out.getValue()
            expected_pack = f"{pack_vers:26.25} - {pack_meta['short_desc']}"
            self.assertEqual(get_pack, expected_pack)
        except:
            macsydata.RemoteModelIndex.list_packages = list_pack
            macsydata.RemoteModelIndex.list_package_vers = list_pack_vers
            macsydata.RemoteModelIndex.get_metadata = meta


    def test_download(self):
        pass

    def test_install(self):
        pass

    def test_uninstall(self):
        pass

    def test_search(self):
        pass

    def test_info(self):
        pass

    def test_list(self):
        pass

    def test_freeze(self):
        pass

    def test_cite(self):
        pass

    def test_check(self):
        pass
