#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
import urllib
import json
import shutil
import tarfile
import glob
import sys

import logging
_log = logging.getLogger(__name__)


class Remote:

    def __init__(self, org="macsy-models"):
        self.org_name = org
        self.base_url = "https://api.github.com"
        self.cache = os.path.join(tempfile.gettempdir(), 'tmp-macsy-cache')

    def _url_json(self, url):
        req = urllib.request.Request(url)
        r = urllib.request.urlopen(req).read()
        j = json.loads(r.decode('utf-8'))
        return j

    def list_packages(self):
        url = f"{self.base_url}/orgs/{self.org_name}/repos"
        _log.debug(f"get {url}")
        packages = self._url_json(url)
        return [p['name'] for p in packages]


    def list_package_vers(self, pack_name):
        url = f"{self.base_url}/repos/{self.org_name}/{pack_name}/tags"
        _log.debug(f"get {url}")
        tags = self._url_json(url)
        return [v['name'] for v in tags]


    def package_download(self, pack_name, vers):
        url = f"{self.base_url}/repos/{self.org_name}/{pack_name}/tarball/{vers}"
        package_cache = os.path.join(self.cache, self.org_name)
        if not os.path.exists(package_cache):
            os.makedirs(package_cache)
        elif os.path.isfile(package_cache):
            raise RuntimeError(f"The tmp cache {package_cache} exist and is a file")
        tmp_archive_path = os.path.join(package_cache, f"{pack_name}-{vers}.tar.gz")

        with urllib.request.urlopen(url) as response, open(tmp_archive_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        return tmp_archive_path


    def unarchive_package(self, path):
        base = os.path.dirname(path)
        *name, vers = '.'.join(os.path.basename(path).split('.')[:-2]).split('-')
        name = '-'.join(name)
        dest_dir = os.path.join(base, name, vers)
        tar = tarfile.open(path, 'r|gz')
        tar.extractall(path=dest_dir)
        src = glob.glob(os.path.join(dest_dir, f"{self.org_name}-{name}-*"))
        if len(src) == 1:
            src = src[0]
        elif len(src) > 1:
            raise RuntimeError(f"Too many matching packages. May be you have to clean {dest_dir}")
        else:
            raise RuntimeError("An error occurred during archive extraction")
        dest = os.path.join(dest_dir, name)
        if os.path.exists(dest):
            shutil.rmtree(dest)
        os.rename(src, dest)
        return dest


class Package:

    def __init__(self, path):
        self.path = path
        self.metadata = self.path.join(self.path, 'metadata')

    def check(path):
        path = os.path.realpath(path)
        if not os.path.exists(path):
            raise RuntimeError()
        elif not os.path.isdir(path):
            raise RuntimeError()
        elif not os.path.exists(os.path.join(path, 'metadata.yml')):
            raise RuntimeError()
        elif not os.path.exists(os.path.join(path, 'definitions')):
            raise RuntimeError()
        elif not os.path.isdir(os.path.join(path, 'definitions')):
            raise RuntimeError()
        elif not os.path.exists(os.path.join(path, 'profiles')):
            raise RuntimeError()
        elif not os.path.isdir(os.path.join(path, 'profiles')):
            raise RuntimeError()
        elif not os.path.exists(os.path.join(path, 'LICENCE')):
            _log.warning("The package {} have not LICENCE file. May be you have not right to use it.")

    # faire la liste de tous les models
    # instancier tous les models


    def help(self):
        readme_path = os.path.join(self.path, 'README')
        with open(readme_path, 'w') as readme:
            for line in readme:
                print(line, file=sys.stderr)

    def _load_metadata(self):
        with open(self.metadata) as raw_metadata:
            metadata = json.loads(raw_metadata.decode('utf-8'))
        return metadata

    def info(self):
        info = """{pack_name} {vers}
author: {auth_name} <{auth_mail}>
{desc}
how to cite:
{}
documentation
{}      
This data are realeased under {}
copyright: {}
"""
        return info

