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
import sys

from macsypy.registries import scan_models_dir, ModelRegistry

from tests import MacsyTest
from macsypy.scripts import macsydata


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
        # some function in macsydata script suppress the traceback
        # but without traceback it's hard to debug test :-(
        sys.tracebacklimit = 1000  # the default value


    def create_fake_package(self, model, definitions=True, profiles=True, metadata=True, readme=True, licence=True):
        pack_path = os.path.join(self.tmpdir, model)
        os.mkdir(pack_path)
        if definitions:
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

        if profiles:
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
        macsydata.RemoteModelIndex.list_package_vers = lambda x, pack: [pack_vers]
        macsydata.RemoteModelIndex.get_metadata = lambda x, pack, vers: pack_meta
        self.create_fake_package('fake_model')
        try:
            with self.catch_io(out=True):
                macsydata.do_available(self.args)
                get_pack = sys.stdout.getvalue().strip()
            pack_name_vers = f"{pack_name} ({pack_vers})"
            expected_pack = f"{pack_name_vers:26.25} - {pack_meta['short_desc']}"
            self.assertEqual(get_pack, expected_pack)
        finally:
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
        pack_name = "nimportnaoik"
        self.args.package = pack_name
        with self.catch_log(log_name='macsydata') as log:
            with self.assertRaises(ValueError):
                macsydata.do_info(self.args)
            log_msg = log.get_value()
        self.assertEqual(log_msg.strip(), f"Models '{pack_name}' not found locally.")

        pack_name = "fake_pack"
        self.args.package = pack_name
        fake_pack_path = self.create_fake_package(pack_name)

        find_local_package = macsydata._find_local_package
        macsydata._find_local_package = lambda x: macsydata.Package(fake_pack_path)
        try:
            with self.catch_io(out=True):
                macsydata.do_info(self.args)
                msg = sys.stdout.getvalue().strip()
        finally:
            macsydata._find_local_package = find_local_package

        expected_info = """fake_pack (0.0b2)

author: auth_name <auth_name@mondomain.fr>

this is a short description of the repos

how to cite:
\t- bla bla
\t- link to publication
\t- ligne 1
\t  ligne 2
\t  ligne 3 et bbbbb
\t  
documentation
\thttp://link/to/the/documentation

This data are released under CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
copyright: 2019, Institut Pasteur, CNRS"""
        self.assertEqual(expected_info, msg)


    def test_list(self):
        fake_packs = ('fake_1', 'fake_2')
        for name in fake_packs:
            self.create_fake_package(name)
        model_dir = self.tmpdir
        registry = ModelRegistry()
        for model_loc in scan_models_dir(model_dir):
            registry.add(model_loc)
        find_all_packages = macsydata._find_all_packages
        macsydata._find_all_packages = lambda: registry

        self.args.verbose = 1
        self.args.outdated = False
        self.args.uptodate = False
        try:
            with self.catch_io(out=True):
                macsydata.do_list(self.args)
                packs = sys.stdout.getvalue().strip()
        finally:
            macsydata._find_all_packages = find_all_packages
        self.assertEqual(packs,
                         "fake_1-0.0b2\nfake_2-0.0b2")


    def test_list_outdated(self):
        fake_packs = ('fake_1', 'fake_2')
        for name in fake_packs:
            self.create_fake_package(name)
        model_dir = self.tmpdir
        registry = ModelRegistry()
        for model_loc in scan_models_dir(model_dir):
            registry.add(model_loc)

        find_all_packages = macsydata._find_all_packages
        macsydata._find_all_packages = lambda: registry

        remote_list_packages_vers = macsydata.RemoteModelIndex.list_package_vers
        macsydata.RemoteModelIndex.list_package_vers = lambda x, name: {'fake_1': ['1.0'],
                                                                        'fake_2': ['0.0b2']}[name]
        self.args.verbose = 1
        self.args.outdated = True
        self.args.uptodate = False

        try:
            with self.catch_io(out=True):
                macsydata.do_list(self.args)
                packs = sys.stdout.getvalue().strip()
        finally:
            macsydata._find_all_packages = find_all_packages
            macsydata.RemoteModelIndex.list_package_vers = remote_list_packages_vers
        self.assertEqual(packs, 'fake_1-1.0 [0.0b2]')


    def test_list_uptodate(self):
        fake_packs = ('fake_1', 'fake_2')
        for name in fake_packs:
            self.create_fake_package(name)
        model_dir = self.tmpdir
        registry = ModelRegistry()
        for model_loc in scan_models_dir(model_dir):
            registry.add(model_loc)
        find_all_packages = macsydata._find_all_packages
        macsydata._find_all_packages = lambda: registry
        remote_list_packages_vers = macsydata.RemoteModelIndex.list_package_vers
        macsydata.RemoteModelIndex.list_package_vers = lambda x, name: {'fake_1': ['1.0'],
                                                                        'fake_2': ['0.0b2']}[name]
        self.args.verbose = 1
        self.args.outdated = False
        self.args.uptodate = True

        try:
            with self.catch_io(out=True):
                macsydata.do_list(self.args)
                packs = sys.stdout.getvalue().strip()
        finally:
            macsydata._find_all_packages = find_all_packages
            macsydata.RemoteModelIndex.list_package_vers = remote_list_packages_vers
        self.assertEqual(packs, 'fake_2-0.0b2')


    def test_list_verbose(self):
        fake_packs = ('fake_1', 'fake_2')
        for name in fake_packs:
            self.create_fake_package(name)
        model_dir = self.tmpdir
        registry = ModelRegistry()
        for model_loc in scan_models_dir(model_dir):
            registry.add(model_loc)
        find_all_packages = macsydata._find_all_packages
        macsydata._find_all_packages = lambda: registry
        remote_list_packages_vers = macsydata.RemoteModelIndex.list_package_vers
        macsydata.RemoteModelIndex.list_package_vers = lambda x, name: {'fake_1': ['1.0'],
                                                                        'fake_2': ['0.0b2']}[name]
        os.unlink(os.path.join(model_dir, 'fake_1', 'metadata.yml'))
        self.args.verbose = 2
        self.args.outdated = False
        self.args.uptodate = False

        try:
            with self.catch_io(out=True):
                with self.catch_log(log_name='macsydata') as log:
                    macsydata.do_list(self.args)
                    log_msg = log.get_value().strip()
                packs = sys.stdout.getvalue().strip()
        finally:
            macsydata._find_all_packages = find_all_packages
            macsydata.RemoteModelIndex.list_package_vers = remote_list_packages_vers
        self.assertEqual(packs, 'fake_2-0.0b2')
        self.assertEqual(log_msg, f"[Errno 2] No such file or directory: '{model_dir}/fake_1/metadata.yml'")


    def test_freeze(self):
        fake_packs = ('fake_1', 'fake_2')
        for name in fake_packs:
            self.create_fake_package(name)
        model_dir = self.tmpdir
        registry = ModelRegistry()
        for model_loc in scan_models_dir(model_dir):
            registry.add(model_loc)
        find_all_packages = macsydata._find_all_packages
        macsydata._find_all_packages = lambda: registry
        try:
            with self.catch_io(out=True):
                macsydata.do_freeze(self.args)
                packs = sys.stdout.getvalue().strip()
        finally:
            macsydata._find_all_packages = find_all_packages
        self.assertEqual(packs,
                         "fake_1==0.0b2\nfake_2==0.0b2")


    def test_cite(self):
        pack_name = "nimportnaoik"
        self.args.package = pack_name
        with self.catch_log(log_name='macsydata') as log:
            with self.assertRaises(ValueError):
                macsydata.do_info(self.args)
            log_msg = log.get_value()
        self.assertEqual(log_msg.strip(), f"Models '{pack_name}' not found locally.")

        pack_name = "fake_pack"
        self.args.package = pack_name
        fake_pack_path = self.create_fake_package(pack_name)

        find_local_package = macsydata._find_local_package
        macsydata._find_local_package = lambda x: macsydata.Package(fake_pack_path)
        try:
            with self.catch_io(out=True):
                macsydata.do_cite(self.args)
                citation = sys.stdout.getvalue().strip()
        finally:
            macsydata._find_local_package = find_local_package
        expected_citation = """To cite fake_pack:

_ bla bla
- link to publication
- ligne 1
  ligne 2
  ligne 3 et bbbbb
  

To cite MacSyFinder:

- Abby SS, Néron B, Ménager H, Touchon M, Rocha EPC (2014)
  MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems.
  PLoS ONE 9(10): e110726. doi:10.1371/journal.pone.0110726"""

        self.assertEqual(expected_citation, citation)


    def test_check(self):
        pack_name = 'fake_1'
        path = self.create_fake_package(pack_name)
        self.args.path = path
        with self.catch_log(log_name='macsydata') as log:
            macsydata.do_check(self.args)
            log_msg = log.get_value().strip()
        expected_msg = f"""If everyone were like you, I'd be out of business
To push the models in organization:
\tcd {os.path.join(self.tmpdir, pack_name)}
Transform the models into a git repository 
\tgit init .
\tgit add .
\tgit commit -m 'initial commit'
add a remote repository to host the models
for instance if you want to add the models to 'macsy-models'
\tgit remote add origin https://github.com/macsy-models/
\tgit tag 0.0b2
\tgit push --tags"""
        self.assertEqual(expected_msg, log_msg)


    def test_check_with_warnings(self):
        pack_name = 'fake_1'
        path = self.create_fake_package(pack_name, readme=False, licence=False)
        self.args.path = path
        with self.catch_log(log_name='macsydata') as log:
            macsydata.do_check(self.args)
            log_msg = log.get_value().strip()
        expected_msg = """The package 'fake_1' have not any LICENCE file. May be you have not right to use it.
The package 'fake_1' have not any README file.
It is better, if you fix warnings above, before to publish these models."""
        self.assertEqual(expected_msg, log_msg)


    def test_check_with_errors(self):
        pack_name = 'fake_1'
        path = self.create_fake_package(pack_name, profiles=False, definitions=False)
        self.args.path = path

        with self.catch_log(log_name='macsydata') as log:
            with self.assertRaises(ValueError):
                macsydata.do_check(self.args)
            log_msg = log.get_value().strip()
        expected_msg = """The package 'fake_1' have no 'definitions' directory.
The package 'fake_1' have no 'profiles' directory.
Please fix issues above, before publishing these models."""
        self.assertEqual(expected_msg, log_msg)


    def test_uninstall(self):
        pack_name = 'fake_1'
        path = self.create_fake_package(pack_name)
        self.args.package = pack_name
        model_loc = scan_models_dir(self.tmpdir)[0]
        find_local_package = macsydata._find_local_package
        macsydata._find_local_package = lambda x: model_loc
        try:
            with self.catch_log(log_name='macsydata') as log:
                macsydata.do_uninstall(self.args)
                log_msg = log.get_value().strip()
        finally:
            macsydata._find_local_package = find_local_package

        expected_msg = f"models '{pack_name}' in {path} uninstalled."
        self.assertEqual(log_msg, expected_msg)
        self.assertFalse(os.path.exists(path))

        self.args.package = 'foo'
        find_local_package = macsydata._find_local_package
        macsydata._find_local_package = lambda x: None
        try:
            with self.catch_log(log_name='macsydata') as log:
                with self.assertRaises(ValueError):
                    macsydata.do_uninstall(self.args)
                log_msg = log.get_value().strip()
        finally:
            macsydata._find_local_package = find_local_package
        expected_msg = f"Models '{self.args.package}' not found locally."
        self.assertEqual(log_msg, expected_msg)
