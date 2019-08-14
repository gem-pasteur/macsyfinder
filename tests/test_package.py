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

import os
import tempfile
import urllib.request
import urllib.error
import json
import io
import shutil
import tarfile
import glob
import yaml
from unittest.mock import patch

from macsypy import package
from macsypy.error import MacsydataError

from tests import MacsyTest


class TestRemote(MacsyTest):

    def setUp(self) -> None:
        self.tmpdir = os.path.join(tempfile.gettempdir(), 'macsy_test_package')
        if os.path.exists(self.tmpdir) and os.path.isdir(self.tmpdir):
            shutil.rmtree(self.tmpdir)
        os.makedirs(self.tmpdir)

    def tearDown(self) -> None:
        try:
            shutil.rmtree(self.tmpdir)
        except:
            pass

    def mocked_requests_get(url):
        class MockResponse:
            def __init__(self, data, status_code):
                self.data = io.BytesIO(bytes(data.encode("utf-8")))
                self.status_code = status_code

            def read(self, length=-1):
                return self.data.read(length)

            def __enter__(self):
                return self

            def __exit__(self, type, value, traceback):
                return False

        if url == 'https://test_url_json/':
            resp = {'fake': ['json', 'response']}
            return MockResponse(json.dumps(resp), 200)
        elif url == 'https://api.github.com/orgs/remote_exists_true':
            resp = {'type': 'Organization'}
            return MockResponse(json.dumps(resp), 200)
        elif url == 'https://api.github.com/orgs/remote_exists_false':
            raise urllib.error.HTTPError(url, 404, 'not found', None, None)
        elif url == 'https://api.github.com/orgs/remote_exists_server_error':
            raise urllib.error.HTTPError(url, 500, 'Server Error', None, None)
        elif url == 'https://api.github.com/orgs/remote_exists_unexpected_error':
            raise urllib.error.HTTPError(url, 204, 'No Content', None, None)
        elif url == 'https://api.github.com/orgs/list_packages/repos':
            resp = [{'name': 'model_1'}, {'name': 'model_2'}]
            return MockResponse(json.dumps(resp), 200)
        elif url == 'https://api.github.com/repos/list_package_vers/model_1/tags':
            resp = [{'name': 'v_1'}, {'name': 'v_2'}]
            return MockResponse(json.dumps(resp), 200)
        elif url == 'https://api.github.com/repos/list_package_vers/model_2/tags':
            raise urllib.error.HTTPError(url, 404, 'not found', None, None)
        elif url == 'https://api.github.com/repos/list_package_vers/model_3/tags':
            raise urllib.error.HTTPError(url, 500, 'Server Error', None, None)
        elif 'https://api.github.com/repos/package_download/fake/tarball/1.0' in url:
            return MockResponse('fake data ' * 2, 200)
        elif url == 'https://api.github.com/repos/package_download/bad_pack/tarball/0.2':
            raise urllib.error.HTTPError(url, 404, 'not found', None, None)
        elif url == 'https://raw.githubusercontent.com/get_metadata/foo/0.0/metadata.yml':
            data = yaml.dump({"author": {"name": "moi"}})
            return MockResponse(data, 200)
        else:
            raise RuntimeError("test non prevu", url)


    def test_init(self):
        rem_exists = package.RemoteModelIndex.remote_exists
        package.RemoteModelIndex.remote_exists = lambda x: True
        try:
            remote = package.RemoteModelIndex()
            remote.cache = self.tmpdir
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists
        self.assertEqual(remote.org_name, 'macsy-models')
        self.assertEqual(remote.base_url, 'https://api.github.com')
        self.assertEqual(remote.cache, self.tmpdir)

        package.RemoteModelIndex.remote_exists = lambda x: True
        try:
            remote = package.RemoteModelIndex(org='foo')
            remote.cache = self.tmpdir
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists
        self.assertEqual(remote.org_name, 'foo')

        package.RemoteModelIndex.remote_exists = lambda x: False
        try:
            with self.assertRaises(ValueError) as ctx:
                package.RemoteModelIndex(org='foo')
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists
        self.assertEqual(str(ctx.exception), "the 'foo' organization does not exist.")


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_url_json(self, mock_urlopen):
        rem_exists = package.RemoteModelIndex.remote_exists
        package.RemoteModelIndex.remote_exists = lambda x: True
        remote = package.RemoteModelIndex(org="nimportnaoik")
        remote.cache = self.tmpdir
        try:
            j = remote._url_json("https://test_url_json/")
            self.assertDictEqual(j, {'fake': ['json', 'response']})
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_remote_exists(self, mock_urlopen):
        remote = package.RemoteModelIndex(org="remote_exists_true")
        remote.cache = self.tmpdir
        exists = remote.remote_exists()
        self.assertTrue(exists)

        remote.org_name = "remote_exists_false"
        exists = remote.remote_exists()
        self.assertFalse(exists)

        remote.org_name = "remote_exists_server_error"
        with self.assertRaises(urllib.error.HTTPError) as ctx:
            remote.remote_exists()
        self.assertEqual(str(ctx.exception),
                         "HTTP Error 500: Server Error")

        remote.org_name = "remote_exists_unexpected_error"
        with self.assertRaises(urllib.error.HTTPError) as ctx:
            remote.remote_exists()
        self.assertEqual(str(ctx.exception),
                         "HTTP Error 204: No Content")

    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_get_metadata(self, mock_urlopen):
        rem_exists = package.RemoteModelIndex.remote_exists
        list_package_vers = package.RemoteModelIndex.list_package_vers
        try:
            vers = '0.0'
            pack_name = 'foo'
            package.RemoteModelIndex.remote_exists = lambda x: True
            package.RemoteModelIndex.list_package_vers = lambda x, pack_name: [vers]
            remote = package.RemoteModelIndex(org="get_metadata")
            remote.cache = self.tmpdir
            metadata = remote.get_metadata(pack_name)
            self.assertDictEqual(metadata, {"author": {"name": "moi"}})
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists
            package.RemoteModelIndex.list_package_vers = list_package_vers

        #################################################
        # The remote package is not versioned (tagged)  #
        #################################################
        try:
            package.RemoteModelIndex.remote_exists = lambda x: True
            package.RemoteModelIndex.list_package_vers = lambda x, pack_name: []
            remote = package.RemoteModelIndex(org="get_metadata")
            with self.assertRaises(MacsydataError) as ctx:
                remote.get_metadata(pack_name)
            self.assertEqual(str(ctx.exception),
                             "No official version available for model 'foo'")
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists
            package.RemoteModelIndex.list_package_vers = list_package_vers

        #####################################
        # The pack version is not available #
        #####################################
        try:
            package.RemoteModelIndex.remote_exists = lambda x: True
            package.RemoteModelIndex.list_package_vers = lambda x, pack_name: ["12"]
            remote = package.RemoteModelIndex(org="get_metadata")
            with self.assertRaises(RuntimeError) as ctx:
                remote.get_metadata(pack_name, vers="1.1")
            self.assertEqual(str(ctx.exception),
                             "The version '1.1' does not exists for model foo.")
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists
            package.RemoteModelIndex.list_package_vers = list_package_vers


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_list_packages(self, mock_urlopen):
        rem_exists = package.RemoteModelIndex.remote_exists
        try:
            package.RemoteModelIndex.remote_exists = lambda x: True
            remote = package.RemoteModelIndex(org="list_packages")
            remote.cache = self.tmpdir
            self.assertListEqual(remote.list_packages(), ['model_1', 'model_2'])
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_list_package_vers(self, mock_urlopen):
        rem_exists = package.RemoteModelIndex.remote_exists
        try:
            package.RemoteModelIndex.remote_exists = lambda x: True
            remote = package.RemoteModelIndex(org="list_package_vers")
            remote.cache = self.tmpdir
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists

        self.assertListEqual(remote.list_package_vers('model_1'), ['v_1', 'v_2'])

        with self.assertRaises(ValueError) as ctx:
            _ = remote.list_package_vers('model_2')
        self.assertEqual(str(ctx.exception), "package 'model_2' does not exists on repos 'list_package_vers'")

        with self.assertRaises(urllib.error.HTTPError) as ctx:
            _ = remote.list_package_vers('model_3')
        self.assertEqual(str(ctx.exception), "HTTP Error 500: Server Error")


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_download(self, mock_urlopen):
        rem_exists = package.RemoteModelIndex.remote_exists
        try:
            package.RemoteModelIndex.remote_exists = lambda x: True
            remote = package.RemoteModelIndex(org="package_download")
            remote.cache = self.tmpdir
            pack_name = "fake"
            pack_vers = "1.0"
            # ensure that remote.cache does not exists
            if os.path.exists(remote.cache):
                if os.path.isdir(remote.cache):
                    shutil.rmtree(remote.cache)
                elif os.path.isfile(remote.cache) or os.path.islink(remote.cache):
                    os.unlink(remote.cache)

            arch_path = remote.download(pack_name, pack_vers)
            self.assertEqual(os.path.join(remote.cache, remote.org_name, f"{pack_name}-{pack_vers}.tar.gz"),
                             arch_path)
            self.assertFileEqual(arch_path, io.StringIO('fake data ' * 2))

            # download again with existing remote.cache and replace old archive
            os.unlink(arch_path)
            arch_path = remote.download(pack_name, pack_vers)
            self.assertFileEqual(arch_path, io.StringIO('fake data ' * 2))

            # download again with existing remote.cache and replace old archive
            os.unlink(arch_path)
            dest = os.path.join(self.tmpdir, 'dest')
            os.makedirs(dest)
            arch_path = remote.download(pack_name, pack_vers, dest=dest)
            self.assertEqual(os.path.join(dest, f'{pack_name}-{pack_vers}.tar.gz'), arch_path)

            # remote cache exist and is a file
            shutil.rmtree(remote.cache)
            open(remote.cache, 'w').close()
            try:
                with self.assertRaises(NotADirectoryError) as ctx:
                    remote.download(pack_name, pack_vers)
                self.assertEqual(str(ctx.exception),
                                 f"The tmp cache '{remote.cache}' already exists")
            finally:
                os.unlink(remote.cache)

            with self.assertRaises(ValueError) as ctx:
                _ = remote.download("bad_pack", "0.2")
            self.assertEqual(str(ctx.exception),
                             "package 'bad_pack-0.2' does not exists on repos 'package_download'")

        finally:
            package.RemoteModelIndex.remote_exists = rem_exists


    def test_unarchive(self):

        def create_pack(dir_, repo, name, vers, key):
            pack_name = f"{name}-{vers}"
            tar_path = os.path.join(dir_, f"{pack_name}.tar.gz")
            with tarfile.open(tar_path, "w:gz") as tar:
                tmp_pack = os.path.join(dir_, f"{repo}-{name}-{key}")
                os.mkdir(tmp_pack)
                for i in range(3):
                    name = f"file_{i}"
                    tmp_file = os.path.join(tmp_pack, name)
                    with open(tmp_file, 'w') as f:
                        f.write(f"Content of file {i}\n")
                tar.add(tmp_pack, arcname=os.path.basename(tmp_pack))
            shutil.rmtree(tmp_pack)
            return tar_path

        pack_name = 'model-toto'
        pack_vers = '2.0'

        rem_exists = package.RemoteModelIndex.remote_exists
        package.RemoteModelIndex.remote_exists = lambda x: True
        try:
            remote = package.RemoteModelIndex(org="package_unarchive")
            arch = create_pack(self.tmpdir, remote.org_name, pack_name, pack_vers, 'e020300')
            remote.cache = self.tmpdir

            model_path = remote.unarchive_package(arch)
            unpacked_path = os.path.join(self.tmpdir, remote.org_name, pack_name, pack_vers, pack_name)
            self.assertEqual(model_path, unpacked_path)
            self.assertTrue(os.path.exists(unpacked_path))
            self.assertTrue(os.path.isdir(unpacked_path))
            self.assertListEqual(sorted(glob.glob(f"{unpacked_path}/*")),
                                 sorted([os.path.join(unpacked_path, f"file_{i}") for i in range(3)])
                                 )
            # create a new archive with old archive not removed
            pack_collision = os.path.join(self.tmpdir, remote.org_name, pack_name, pack_vers,
                                          f"{remote.org_name}-{pack_name}-e022222")
            os.mkdir(pack_collision)
            arch = create_pack(self.tmpdir, remote.org_name, pack_name, pack_vers, 'e0203333')
            with self.assertRaises(MacsydataError) as ctx:
                remote.unarchive_package(arch)
            self.assertEqual(str(ctx.exception),
                             f"Too many matching packages. May be you have to clean "
                             f"{os.path.join(self.tmpdir, remote.org_name, pack_name, pack_vers)}")
        finally:
            package.RemoteModelIndex.remote_exists = rem_exists


class TestPackage(MacsyTest):

    def setUp(self) -> None:
        self.tmpdir = os.path.join(tempfile.gettempdir(), 'macsy_test_package')
        if os.path.exists(self.tmpdir) and os.path.isdir(self.tmpdir):
            shutil.rmtree(self.tmpdir)
        os.makedirs(self.tmpdir)

        self.metadata = {"author": {"name": "auth_name",
                                    "email": "auth_name@mondomain.fr"},
                         "short_desc": "this is a short description of the repos",
                         "vers": "0.0b2",
                         "cite": ["bla bla",
                                  "link to publication",
                                  """ligne 1
ligne 2
ligne 3 et bbbbb
"""],
                         "doc": "http://link/to/the/documentation",
                         "licence": "CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)",
                         "copyright": "2019, Institut Pasteur, CNRS"
                          }


    def tearDown(self) -> None:
        try:
            shutil.rmtree(self.tmpdir)
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


    def test_init(self):
        fake_pack_path = self.create_fake_package('fake_model')
        pack = package.Package(fake_pack_path)
        self.assertEqual(pack.path, fake_pack_path)
        self.assertEqual(pack.readme, os.path.join(fake_pack_path, 'README'))
        self.assertEqual(pack.name, 'fake_model')
        self.assertEqual(pack.metadata_path, os.path.join(fake_pack_path, 'metadata.yml'))


    def test_metadata(self):
        fake_pack_path = self.create_fake_package('fake_model')
        pack = package.Package(fake_pack_path)
        self.assertDictEqual(pack.metadata, self.metadata)
        self.assertDictEqual(pack.metadata, self.metadata)

    def test_find_readme(self):
        fake_pack_path = self.create_fake_package('fake_model')
        pack = package.Package(fake_pack_path)
        for ext in ('', '.rst', '.md'):
            readme_path = os.path.join(pack.path, 'README' + ext)
            os.rename(pack.readme, readme_path)
            pack.readme = readme_path
            self.assertEqual(pack._find_readme(), readme_path)
        readme_path = os.path.join(pack.path, 'README.foo')
        os.rename(pack.readme, readme_path)
        self.assertIsNone(pack._find_readme())


    def test_check_structure(self):
        fake_pack_path = self.create_fake_package('fake_model')
        pack = package.Package(fake_pack_path)
        errors, warnings = pack._check_structure()
        self.assertListEqual(errors, [])
        self.assertListEqual(warnings, [])


    def test_check_structure_bad_path(self):
        foobar = os.path.join(self.tmpdir, "foobar")
        pack = package.Package(foobar)
        errors, warnings = pack._check_structure()
        self.assertListEqual(errors, ["The package 'foobar' does not exists."])
        self.assertListEqual(warnings, [])

        open(foobar, 'w').close()
        errors, warnings = pack._check_structure()
        self.assertListEqual(errors, ["'foobar' is not a directory "])
        self.assertListEqual(warnings, [])


    def test_check_structure_no_def(self):
        fake_pack_path = self.create_fake_package('fake_model', xml=False)
        pack = package.Package(fake_pack_path)
        errors, warnings = pack._check_structure()

        self.assertListEqual(errors, ["The package 'fake_model' have no 'definitions' directory."])
        self.assertListEqual(warnings, [])

        open(os.path.join(pack.path, 'definitions'), 'w').close()
        errors, warnings = pack._check_structure()
        self.assertListEqual(errors, ["'/tmp/macsy_test_package/fake_model/definitions' is not a directory."])
        self.assertListEqual(warnings, [])


    def test_check_structure_no_profiles(self):
        fake_pack_path = self.create_fake_package('fake_model', hmm=False)
        pack = package.Package(fake_pack_path)
        errors, warnings = pack._check_structure()

        self.assertListEqual(errors, ["The package 'fake_model' have no 'profiles' directory."])
        self.assertListEqual(warnings, [])

        open(os.path.join(pack.path, 'profiles'), 'w').close()
        errors, warnings = pack._check_structure()
        self.assertListEqual(errors, ["'/tmp/macsy_test_package/fake_model/profiles' is not a directory."])
        self.assertListEqual(warnings, [])


    def test_check_structure_no_metadata(self):
        fake_pack_path = self.create_fake_package('fake_model', metadata=False)
        pack = package.Package(fake_pack_path)
        errors, warnings = pack._check_structure()

        self.assertListEqual(errors, ["The package 'fake_model' have no 'metadata.yml'."])
        self.assertListEqual(warnings, [])


    def test_check_structure_no_readme(self):
        fake_pack_path = self.create_fake_package('fake_model', readme=False)
        pack = package.Package(fake_pack_path)
        errors, warnings = pack._check_structure()

        self.assertEqual(errors, [])
        self.assertEqual(warnings, ["The package 'fake_model' have not any README file."])


    def test_check_structure_no_licence(self):
        fake_pack_path = self.create_fake_package('fake_model', licence=False)
        pack = package.Package(fake_pack_path)
        errors, warnings = pack._check_structure()

        self.assertEqual(errors, [])
        self.assertEqual(warnings, ["The package 'fake_model' have not any LICENCE file. "
                                    "May be you have not right to use it."])


    def test_check_metadata(self):
        fake_pack_path = self.create_fake_package('fake_model')
        pack = package.Package(fake_pack_path)
        errors, warnings = pack._check_metadata()
        self.assertEqual(errors, [])
        self.assertEqual(warnings, [])

        load_metadata_meth = package.Package._load_metadata

        #############
        # No author #
        #############
        no_auth_meta_data = self.metadata.copy()
        del no_auth_meta_data['author']
        try:
            package.Package._load_metadata = lambda x: no_auth_meta_data
            pack = package.Package(fake_pack_path)
            errors, warnings = pack._check_metadata()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertEqual(errors, ["field 'author' is mandatory in metadata_path."])
        self.assertEqual(warnings, [])

        #################
        # No short desc #
        #################
        no_short_desc_metadata = self.metadata.copy()
        del no_short_desc_metadata['short_desc']
        try:
            package.Package._load_metadata = lambda x: no_short_desc_metadata
            pack = package.Package(fake_pack_path)
            errors, warnings = pack._check_metadata()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertEqual(errors, ["field 'short_desc' is mandatory in metadata_path."])
        self.assertEqual(warnings, [])

        ###########
        # No vers #
        ###########
        no_vers_metadata = self.metadata.copy()
        del no_vers_metadata['vers']
        try:
            package.Package._load_metadata = lambda x: no_vers_metadata
            pack = package.Package(fake_pack_path)
            errors, warnings = pack._check_metadata()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertEqual(errors, ["field 'vers' is mandatory in metadata_path."])
        self.assertEqual(warnings, [])

        ###########
        # No cite #
        ###########
        no_cite_metadata = self.metadata.copy()
        del no_cite_metadata['cite']
        try:
            package.Package._load_metadata = lambda x: no_cite_metadata
            pack = package.Package(fake_pack_path)
            errors, warnings = pack._check_metadata()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertEqual(errors, [])
        self.assertEqual(warnings, ["It's better if the field 'cite' is setup in metadata_path file"])

        ##########
        # No doc #
        ##########
        no_doc_metadata = self.metadata.copy()
        del no_doc_metadata['doc']
        try:
            package.Package._load_metadata = lambda x: no_doc_metadata
            pack = package.Package(fake_pack_path)
            errors, warnings = pack._check_metadata()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertEqual(errors, [])
        self.assertEqual(warnings, ["It's better if the field 'doc' is setup in metadata_path file"])

        ##############
        # No licence #
        ##############
        no_licence_metadata = self.metadata.copy()
        del no_licence_metadata['licence']
        try:
            package.Package._load_metadata = lambda x: no_licence_metadata
            pack = package.Package(fake_pack_path)
            errors, warnings = pack._check_metadata()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertEqual(errors, [])
        self.assertEqual(warnings, ["It's better if the field 'licence' is setup in metadata_path file"])

        ################
        # No copyright #
        ################
        no_copyright_metadata = self.metadata.copy()
        del no_copyright_metadata['copyright']
        try:
            package.Package._load_metadata = lambda x: no_copyright_metadata
            pack = package.Package(fake_pack_path)
            errors, warnings = pack._check_metadata()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertEqual(errors, [])
        self.assertEqual(warnings, ["It's better if the field 'copyright' is setup in metadata_path file"])

        ##################
        # No author name #
        ##################
        # this test must the last of the set
        # because we remove a key in author value
        # the copy is a shallow copy
        # so author value is a reference to the good_metadata[author]
        # side effect
        no_auth_name_meta_data = self.metadata.copy()
        del no_auth_name_meta_data['author']['name']
        try:
            package.Package._load_metadata = lambda x: no_auth_name_meta_data
            pack = package.Package(fake_pack_path)
            errors, warnings = pack._check_metadata()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertEqual(errors, ["field 'author.name' is mandatory in metadata_path."])
        self.assertEqual(warnings, [])


    def test_check(self):
        fake_pack_path = self.create_fake_package('fake_model')
        load_metadata_meth = package.Package._load_metadata
        bad_meta_data = {"short_desc": "this is a short description of the repos",
                         "doc": "http://link/to/the/documentation",
                         "copyright": "2019, Institut Pasteur, CNRS"
                         }
        try:
            package.Package._load_metadata = lambda x: bad_meta_data
            pack = package.Package(fake_pack_path)
            errors, warnings = pack.check()
        finally:
            package.Package._load_metadata = load_metadata_meth
        self.assertListEqual(errors,
                             ["field 'author' is mandatory in metadata_path.",
                              "field 'vers' is mandatory in metadata_path."])
        self.assertListEqual(warnings,
                             ["It's better if the field 'cite' is setup in metadata_path file",
                              "It's better if the field 'licence' is setup in metadata_path file"])

    def test_help(self):
        fake_pack_path = self.create_fake_package('fake_model', licence=False)
        pack = package.Package(fake_pack_path)

        receive_help = pack.help()
        self.assertEqual(receive_help, "# This a README\n")

        os.unlink(os.path.join(fake_pack_path, 'README'))
        pack = package.Package(fake_pack_path)
        receive_help = pack.help()
        self.assertEqual(receive_help, "No help available for package 'fake_model'.")


    def test_info(self):
        fake_pack_path = self.create_fake_package('fake_model', licence=False)
        pack = package.Package(fake_pack_path)

        info = pack.info()
        expected_info = """
fake_model (0.0b2)

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
copyright: 2019, Institut Pasteur, CNRS
"""
        self.assertEqual(info, expected_info)

        load_metadata_meth = package.Package._load_metadata
        ###########
        # No cite #
        ###########
        no_cite_metadata = self.metadata.copy()
        del no_cite_metadata['cite']
        try:
            package.Package._load_metadata = lambda x: no_cite_metadata
            pack = package.Package(fake_pack_path)
            info = pack.info()
        finally:
            package.Package._load_metadata = load_metadata_meth
        expected_info = """
fake_model (0.0b2)

author: auth_name <auth_name@mondomain.fr>

this is a short description of the repos

how to cite:
\t- No citation available
\t  
documentation
\thttp://link/to/the/documentation

This data are released under CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
copyright: 2019, Institut Pasteur, CNRS
"""
        self.assertEqual(info, expected_info)

        ##########
        # No doc #
        ##########
        no_doc_metadata = self.metadata.copy()
        del no_doc_metadata['doc']
        try:
            package.Package._load_metadata = lambda x: no_doc_metadata
            pack = package.Package(fake_pack_path)
            info = pack.info()
        finally:
            package.Package._load_metadata = load_metadata_meth
        expected_info = """
fake_model (0.0b2)

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
\tNo documentation available

This data are released under CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
copyright: 2019, Institut Pasteur, CNRS
"""
        self.assertEqual(info, expected_info)

        ##############
        # No licence #
        ##############
        no_licence_metadata = self.metadata.copy()
        del no_licence_metadata['licence']
        try:
            package.Package._load_metadata = lambda x: no_licence_metadata
            pack = package.Package(fake_pack_path)
            info = pack.info()
        finally:
            package.Package._load_metadata = load_metadata_meth
        expected_info = """
fake_model (0.0b2)

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

This data are released under No licence available
copyright: 2019, Institut Pasteur, CNRS
"""
        self.assertEqual(info, expected_info)
