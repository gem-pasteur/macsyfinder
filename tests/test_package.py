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
from unittest.mock import patch, MagicMock, Mock

from macsypy import package

from tests import MacsyTest


class TestRemote(MacsyTest):

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


    @patch('urllib.request.urlopen')
    def test_url_json(self, mock_urlopen):
        cm = MagicMock()
        cm.getcode.return_value = 200
        resp = {'fake': ['json', 'response']}
        cm.read.return_value = bytes(json.dumps(resp).encode(encoding='utf-8'))
        cm.__enter__.return_value = cm
        mock_urlopen.return_value = cm

        rem_exists = package.Remote.remote_exists
        package.Remote.remote_exists = lambda x: True
        remote = package.Remote(org="nimportnaoik")
        try:
            j = remote._url_json("https://nimportnaoik")
            self.assertDictEqual(j, resp)
        finally:
            package.Remote.remote_exists = rem_exists


    @patch('urllib.request.urlopen')
    def test_remote_exists(self, mock_urlopen):
        cm = MagicMock()
        cm.getcode.return_value = 200
        resp = {'type': 'Organization'}
        cm.read.return_value = bytes(json.dumps(resp).encode(encoding='utf-8'))
        cm.__enter__.return_value = cm
        mock_urlopen.return_value = cm

        remote = package.Remote(org="nimportnaoik")
        exists = remote.remote_exists()
        self.assertTrue(exists)

        def raise_404_http_error(_, url):
            raise urllib.error.HTTPError(url, 404, 'not found', None, None)

        url_json = package.Remote._url_json
        package.Remote._url_json = raise_404_http_error
        try:
            exists = remote.remote_exists()
        finally:
            package.Remote._url_json = url_json
        self.assertFalse(exists)

        def raise_503_http_error(_, url):
            raise urllib.error.HTTPError(url, 503, 'server error', None, None)

        url_json = package.Remote._url_json
        package.Remote._url_json = raise_503_http_error
        try:
            self.assertRaises(urllib.error.HTTPError, remote.remote_exists)
        finally:
            package.Remote._url_json = url_json


    def test_list_packages(self):
        rem_exists = package.Remote.remote_exists
        package.Remote.remote_exists = lambda x: True
        try:
            remote = package.Remote(org="nimportnaoik")
        finally:
            package.Remote.remote_exists = rem_exists

        url_json = package.Remote._url_json
        resp = [{'name': 'model_1'}, {'name': 'model_2'}]
        package.Remote._url_json = lambda x, url: resp
        try:
            self.assertListEqual(remote.list_packages(), ['model_1', 'model_2'])
        finally:
            package.Remote._url_json = url_json


    def test_list_package_vers(self):
        rem_exists = package.Remote.remote_exists
        package.Remote.remote_exists = lambda x: True
        try:
            remote = package.Remote(org="nimportnaoik")
        finally:
            package.Remote.remote_exists = rem_exists

        url_json = package.Remote._url_json
        resp = [{'name': 'v_1'}, {'name': 'v_2'}]
        package.Remote._url_json = lambda x, url: resp
        try:
            self.assertListEqual(remote.list_package_vers('model_1'), ['v_1', 'v_2'])
        finally:
            package.Remote._url_json = url_json

        def raise_404_http_error(_, url):
            raise urllib.error.HTTPError(url, 404, 'not found', None, None)

        url_json = package.Remote._url_json
        package.Remote._url_json = raise_404_http_error
        try:
            with self.assertRaises(RuntimeError) as ctx:
                _ = remote.list_package_vers('model_1')
            self.assertEqual(str(ctx.exception), "package 'model_1' does not exists on repos 'nimportnaoik'")
        finally:
            package.Remote._url_json = url_json


class TestPackage(MacsyTest):
    pass