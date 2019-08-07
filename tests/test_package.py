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
from unittest.mock import patch

from macsypy import package

from tests import MacsyTest


# class TestRemote(MacsyTest):
#
#     def setUp(self) -> None:
#         self.tmpdir = os.path.join(tempfile.gettempdir(), 'macsy_test_package')
#         if os.path.exists(self.tmpdir) and os.path.isdir(self.tmpdir):
#             shutil.rmtree(self.tmpdir)
#         os.makedirs(self.tmpdir)
#
#     def tearDown(self) -> None:
#         try:
#             shutil.rmtree(self.tmpdir)
#         except:
#             pass
#
#     def mocked_requests_get(url):
#         class MockResponse:
#             def __init__(self, data, status_code):
#                 self.data = io.BytesIO(bytes(data.encode("utf-8")))
#                 self.status_code = status_code
#
#             def read(self, length=-1):
#                 return self.data.read(length)
#
#             def __enter__(self):
#                 return self
#
#             def __exit__(self, type, value, traceback):
#                 return False
#         if url == 'https://test_url_json/':
#             resp = {'fake': ['json', 'response']}
#             return MockResponse(json.dumps(resp), 200)
#         elif url == 'https://api.github.com/orgs/remote_exists_true':
#             resp = {'type': 'Organization'}
#             return MockResponse(json.dumps(resp), 200)
#         elif url == 'https://api.github.com/orgs/remote_exists_false':
#             raise urllib.error.HTTPError(url, 404, 'not found', None, None)
#         elif url == 'https://api.github.com/orgs/remote_exists_server_error':
#             raise urllib.error.HTTPError(url, 500, 'Server Error', None, None)
#         elif url == 'https://api.github.com/orgs/remote_exists_unexpected_error':
#             raise urllib.error.HTTPError(url, 204, 'No Content', None, None)
#         elif url == 'https://api.github.com/orgs/list_packages/repos':
#             resp = [{'name': 'model_1'}, {'name': 'model_2'}]
#             return MockResponse(json.dumps(resp), 200)
#         elif url == 'https://api.github.com/repos/list_package_vers/model_1/tags':
#             resp = [{'name': 'v_1'}, {'name': 'v_2'}]
#             return MockResponse(json.dumps(resp), 200)
#         elif url == 'https://api.github.com/repos/list_package_vers/model_2/tags':
#             raise urllib.error.HTTPError(url, 404, 'not found', None, None)
#         elif url == 'https://api.github.com/repos/list_package_vers/model_3/tags':
#             raise urllib.error.HTTPError(url, 500, 'Server Error', None, None)
#         elif 'https://api.github.com/repos/package_download/fake/tarball/1.0' in url:
#             return MockResponse('fake data ' * 2, 200)
#         elif 'https://api.github.com/repos/package_download/bad_pack/tarball/name':
#             raise urllib.error.HTTPError(url, 404, 'not found', None, None)
#         else:
#             raise RuntimeError("test non prevu", args)
#
#
#     def test_init(self):
#         rem_exists = package.Remote.remote_exists
#         package.Remote.remote_exists = lambda x: True
#         try:
#             remote = package.Remote()
#             remote.cache = self.tmpdir
#         finally:
#             package.Remote.remote_exists = rem_exists
#         self.assertEqual(remote.org_name, 'macsy-models')
#         self.assertEqual(remote.base_url, 'https://api.github.com')
#         self.assertEqual(remote.cache, self.tmpdir)
#
#         package.Remote.remote_exists = lambda x: True
#         try:
#             remote = package.Remote(org='foo')
#             remote.cache = self.tmpdir
#         finally:
#             package.Remote.remote_exists = rem_exists
#         self.assertEqual(remote.org_name, 'foo')
#
#
#     @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
#     def test_url_json(self, mock_urlopen):
#         rem_exists = package.Remote.remote_exists
#         package.Remote.remote_exists = lambda x: True
#         remote = package.Remote(org="nimportnaoik")
#         remote.cache = self.tmpdir
#         try:
#             j = remote._url_json("https://test_url_json/")
#             self.assertDictEqual(j, {'fake': ['json', 'response']})
#         finally:
#             package.Remote.remote_exists = rem_exists
#
#
#     @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
#     def test_remote_exists(self, mock_urlopen):
#         remote = package.Remote(org="remote_exists_true")
#         remote.cache = self.tmpdir
#         exists = remote.remote_exists()
#         self.assertTrue(exists)
#
#         remote.org_name = "remote_exists_false"
#         exists = remote.remote_exists()
#         self.assertFalse(exists)
#
#         remote.org_name = "remote_exists_server_error"
#         with self.assertRaises(urllib.error.HTTPError) as ctx:
#             remote.remote_exists()
#         self.assertEqual(str(ctx.exception),
#                          "HTTP Error 500: Server Error")
#
#         remote.org_name = "remote_exists_unexpected_error"
#         with self.assertRaises(urllib.error.HTTPError) as ctx:
#             remote.remote_exists()
#         self.assertEqual(str(ctx.exception),
#                          "HTTP Error 204: No Content")
#
#
#     @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
#     def test_list_packages(self, mock_urlopen):
#         rem_exists = package.Remote.remote_exists
#         try:
#             package.Remote.remote_exists = lambda x: True
#             remote = package.Remote(org="list_packages")
#             remote.cache = self.tmpdir
#             self.assertListEqual(remote.list_packages(), ['model_1', 'model_2'])
#         finally:
#             package.Remote.remote_exists = rem_exists
#
#
#     @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
#     def test_list_package_vers(self, mock_urlopen):
#         rem_exists = package.Remote.remote_exists
#         try:
#             package.Remote.remote_exists = lambda x: True
#             remote = package.Remote(org="list_package_vers")
#             remote.cache = self.tmpdir
#         finally:
#             package.Remote.remote_exists = rem_exists
#
#         self.assertListEqual(remote.list_package_vers('model_1'), ['v_1', 'v_2'])
#
#         with self.assertRaises(RuntimeError) as ctx:
#             _ = remote.list_package_vers('model_2')
#         self.assertEqual(str(ctx.exception), "package 'model_2' does not exists on repos 'list_package_vers'")
#
#         with self.assertRaises(urllib.error.HTTPError) as ctx:
#             _ = remote.list_package_vers('model_3')
#         self.assertEqual(str(ctx.exception), "HTTP Error 500: Server Error")
#
#
#     @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
#     def test_package_download(self, mock_urlopen):
#         rem_exists = package.Remote.remote_exists
#         try:
#             package.Remote.remote_exists = lambda x: True
#             remote = package.Remote(org="package_download")
#             remote.cache = self.tmpdir
#             pack_name = "fake"
#             pack_vers = "1.0"
#             # ensure that remote cache does not exists
#             if os.path.exists(remote.cache):
#                 if os.path.isdir(remote.cache):
#                     shutil.rmtree(remote.cache)
#                 elif os.path.isfile(remote.cache) or os.path.islink(remote.cache):
#                     os.unlink(remote.cache)
#
#             arch_path = remote.package_download(pack_name, pack_vers)
#             self.assertEqual(os.path.join(remote.cache, remote.org_name, f"{pack_name}-{pack_vers}.tar.gz"),
#                              arch_path)
#             self.assertFileEqual(arch_path, io.StringIO('fake data ' * 2))
#
#             # download again with existing remote cache and replace old archive
#             os.unlink(arch_path)
#             arch_path = remote.package_download(pack_name, pack_vers)
#             self.assertFileEqual(arch_path, io.StringIO('fake data ' * 2))
#
#             # remote cache exist and is a file
#             shutil.rmtree(remote.cache)
#             open(remote.cache, 'w').close()
#             try:
#                 with self.assertRaises(NotADirectoryError) as ctx:
#                     remote.package_download(pack_name, pack_vers)
#                 self.assertEqual(str(ctx.exception),
#                                  f"The tmp cache '{remote.cache}' already exists")
#             finally:
#                 os.unlink(remote.cache)
#
#             with self.assertRaises(RuntimeError) as ctx:
#                 _ = remote.package_download("bad_pack", "name")
#             self.assertEqual(str(ctx.exception),
#                              "package 'bad_pack-name' does not exists on repos 'package_download'")
#
#         finally:
#             package.Remote.remote_exists = rem_exists
#
#
#     def test_unarchive(self):
#
#         def create_pack(dir_, repo, name, vers, key):
#             pack_name = f"{name}-{vers}"
#             tar_path = os.path.join(dir_, f"{pack_name}.tar.gz")
#             with tarfile.open(tar_path, "w:gz") as tar:
#                 tmp_pack = os.path.join(dir_, f"{repo}-{name}-{key}")
#                 os.mkdir(tmp_pack)
#                 for i in range(3):
#                     name = f"file_{i}"
#                     tmp_file = os.path.join(tmp_pack, name)
#                     with open(tmp_file, 'w') as f:
#                         f.write(f"Content of file {i}\n")
#                 tar.add(tmp_pack, arcname=os.path.basename(tmp_pack))
#             shutil.rmtree(tmp_pack)
#             return tar_path
#
#         pack_name = 'model-toto'
#         pack_vers = '2.0'
#
#         rem_exists = package.Remote.remote_exists
#         package.Remote.remote_exists = lambda x: True
#         try:
#             remote = package.Remote(org="package_unarchive")
#             arch = create_pack(self.tmpdir, remote.org_name, pack_name, pack_vers, 'e020300')
#             remote.cache = self.tmpdir
#
#             model_path = remote.unarchive_package(arch)
#             unpacked_path = os.path.join(self.tmpdir, remote.org_name, pack_name, pack_vers, pack_name)
#             self.assertEqual(model_path, unpacked_path)
#             self.assertTrue(os.path.exists(unpacked_path))
#             self.assertTrue(os.path.isdir(unpacked_path))
#             self.assertListEqual(sorted(glob.glob(f"{unpacked_path}/*")),
#                                  sorted([os.path.join(unpacked_path, f"file_{i}") for i in range(3)])
#                                  )
#             # create a new archive with old archive not removed
#             pack_collision = os.path.join(self.tmpdir, remote.org_name, pack_name, pack_vers,
#                                           f"{remote.org_name}-{pack_name}-e022222")
#             os.mkdir(pack_collision)
#             arch = create_pack(self.tmpdir, remote.org_name, pack_name, pack_vers, 'e0203333')
#             with self.assertRaises(RuntimeError) as ctx:
#                 remote.unarchive_package(arch)
#             self.assertEqual(str(ctx.exception),
#                              f"Too many matching packages. May be you have to clean "
#                              f"{os.path.join(self.tmpdir, remote.org_name, pack_name, pack_vers)}")
#         finally:
#             package.Remote.remote_exists = rem_exists

class TestPackage(MacsyTest):

    def setUp(self) -> None:
        self.tmpdir = os.path.join(tempfile.gettempdir(), 'macsy_test_package')
        if os.path.exists(self.tmpdir) and os.path.isdir(self.tmpdir):
            shutil.rmtree(self.tmpdir)
        os.makedirs(self.tmpdir)

    def tearDown(self) -> None:
        try:
            #shutil.rmtree(self.tmpdir)
            pass
        except:
            pass

    def create_fake_package(self, model,
                            metadata=True,
                            readme=True,
                            licence=True):
        pack_path = os.path.join(self.tmpdir, model)
        os.mkdir(pack_path)
        for name, ext in (('definitions', 'xml'), ('profiles', 'hmm')):
            sub_dir_path = os.path.join(pack_path, name)
            os.mkdir(sub_dir_path)
            for i in range(2):
                open(os.path.join(sub_dir_path, f"file_{i}.{ext}"), 'w').close()
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
        fake_pack = self.create_fake_package('fake_model')
        pack = package.Package(fake_pack)
        self.assertEqual(pack.path, fake_pack)