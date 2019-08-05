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
from unittest.mock import patch, MagicMock, Mock

from macsypy import package

from tests import MacsyTest


class TestRemote(MacsyTest):

    def mocked_requests_get(*args, **kwargs):
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
        if args[0] == 'https://test_url_json/':
            resp = {'fake': ['json', 'response']}
            return MockResponse(json.dumps(resp), 200)
        elif args[0] == 'https://api.github.com/orgs/remote_exists_true':
            resp = {'type': 'Organization'}
            return MockResponse(json.dumps(resp), 200)
        elif args[0] == 'https://api.github.com/orgs/remote_exists_false':
            raise urllib.error.HTTPError(args[0], 404, 'not found', None, None)
        elif args[0] == 'https://api.github.com/orgs/remote_exists_server_error':
            raise urllib.error.HTTPError(args[0], 503, 'not found', None, None)
        elif args[0] == 'https://api.github.com/orgs/list_packages/repos':
            resp = [{'name': 'model_1'}, {'name': 'model_2'}]
            return MockResponse(json.dumps(resp), 200)
        elif args[0] == 'https://api.github.com/repos/list_package_vers/model_1/tags':
            resp = [{'name': 'v_1'}, {'name': 'v_2'}]
            return MockResponse(json.dumps(resp), 200)
        elif args[0] == 'https://api.github.com/repos/list_package_vers/model_2/tags':
            raise urllib.error.HTTPError(args[0], 404, 'not found', None, None)
        elif 'https://api.github.com/repos/package_download/fake/tarball/1.0' in args[0]:
            return MockResponse('fake data ' * 2, 200)
        elif 'https://api.github.com/repos/package_download/bad_pack/tarball/name':
            raise urllib.error.HTTPError(args[0], 404, 'not found', None, None)
        else:
            raise RuntimeError("test non prevu", args)



    def test_init(self):
        rem_exists = package.Remote.remote_exists
        package.Remote.remote_exists = lambda x: True
        try:
            remote = package.Remote()
        finally:
            package.Remote.remote_exists = rem_exists
        self.assertEqual(remote.org_name, 'macsy-models')
        self.assertEqual(remote.base_url, 'https://api.github.com')
        self.assertEqual(remote.cache, os.path.join(tempfile.gettempdir(), 'tmp-macsy-cache'))

        package.Remote.remote_exists = lambda x: True
        try:
            remote = package.Remote(org='foo')
        finally:
            package.Remote.remote_exists = rem_exists
        self.assertEqual(remote.org_name, 'foo')


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_url_json(self, mock_urlopen):
        rem_exists = package.Remote.remote_exists
        package.Remote.remote_exists = lambda x: True
        remote = package.Remote(org="nimportnaoik")
        try:
            j = remote._url_json("https://test_url_json/")
            self.assertDictEqual(j, {'fake': ['json', 'response']})
        finally:
            package.Remote.remote_exists = rem_exists


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_remote_exists(self, mock_urlopen):
        remote = package.Remote(org="remote_exists_true")
        exists = remote.remote_exists()
        self.assertTrue(exists)

        remote.org_name = "remote_exists_false"
        exists = remote.remote_exists()
        self.assertFalse(exists)

        remote.org_name = "remote_exists_server_error"
        self.assertRaises(urllib.error.HTTPError, remote.remote_exists)


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_list_packages(self, mock_urlopen):
        rem_exists = package.Remote.remote_exists
        try:
            package.Remote.remote_exists = lambda x: True
            remote = package.Remote(org="list_packages")
            self.assertListEqual(remote.list_packages(), ['model_1', 'model_2'])
        finally:
            package.Remote.remote_exists = rem_exists

    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_list_package_vers(self, mock_urlopen):
        rem_exists = package.Remote.remote_exists
        try:
            package.Remote.remote_exists = lambda x: True
            remote = package.Remote(org="list_package_vers")
        finally:
            package.Remote.remote_exists = rem_exists

        self.assertListEqual(remote.list_package_vers('model_1'), ['v_1', 'v_2'])

        with self.assertRaises(RuntimeError) as ctx:
            _ = remote.list_package_vers('model_2')
        self.assertEqual(str(ctx.exception), "package 'model_2' does not exists on repos 'list_package_vers'")


    @patch('urllib.request.urlopen', side_effect=mocked_requests_get)
    def test_package_download(self, mock_urlopen):
        rem_exists = package.Remote.remote_exists
        try:
            package.Remote.remote_exists = lambda x: True
            remote = package.Remote(org="package_download")
            pack_name = "fake"
            pack_vers = "1.0"
            arch_path = remote.package_download(pack_name, pack_vers)
            self.assertEqual(os.path.join(remote.cache, remote.org_name, f"{pack_name}-{pack_vers}.tar.gz"),
                             arch_path)
            self.assertFileEqual(arch_path, io.StringIO('fake data ' * 2))
            with self.assertRaises(RuntimeError) as ctx:
                _ = remote.package_download("bad_pack", "name")
            self.assertEqual(str(ctx.exception),
                             "package 'bad_pack-name' does not exists on repos 'package_download'")
        finally:
            package.Remote.remote_exists = rem_exists




class TestPackage(MacsyTest):
    pass