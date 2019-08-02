#!/usr/bin/env python
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
import json
import yaml
import shutil
import tarfile
import glob
import sys
from typing import List, Dict, Any

import logging
_log = logging.getLogger(__name__)

from .config import NoneConfig
from .registries import ModelLocation
from .definition_parser import DefinitionParser
from .model import ModelBank
from .gene import GeneBank, ProfileFactory


class Remote:

    def __init__(self, org="macsy-models"):
        """

        :param org: The name of the organization on github where are stored the models
        """
        self.org_name = org
        self.base_url = "https://api.github.com"
        self.cache = os.path.join(tempfile.gettempdir(), 'tmp-macsy-cache')
        self.remote_exists()


    def _url_json(self, url: str) -> Dict:
        """
        Get the url, deserialize the data as json

        :param str url: the url to dowload
        :return: the json corresponding to the response url
        """
        req = urllib.request.Request(url)
        r = urllib.request.urlopen(req).read()
        j = json.loads(r.decode('utf-8'))
        return j


    def remote_exists(self) -> bool:
        """
        check if the remote exists and is an organization
        :return:
        """
        try:
            url = f"{self.base_url}/orgs/{self.org_name}"
            _log.debug(f"get {url}")
            remote = self._url_json(url)
            return remote["type"] == 'Organization'
        except urllib.error.HTTPError as err:
            if 400 <= err.code < 500:
                return False
            elif err.code >= 500:
                raise err from None
            else:
                raise err from None


    def list_packages(self) -> List[str]:
        """
        list all model packages availables on a model repos
        :return: The list of package names.
        """
        url = f"{self.base_url}/orgs/{self.org_name}/repos"
        _log.debug(f"get {url}")
        packages = self._url_json(url)
        return [p['name'] for p in packages]


    def list_package_vers(self, pack_name: str) -> List:
        """
        List all available versions from github model repos for a given package

        :param str pack_name: the name of the package
        :return: the list of the versions
        """
        url = f"{self.base_url}/repos/{self.org_name}/{pack_name}/tags"
        _log.debug(f"get {url}")
        try:
            tags = self._url_json(url)
        except urllib.error.HTTPError as err:
            if 400 <= err.code < 500:
                raise RuntimeError(f"package '{pack_name}'' does not exists on repos '{self.org_name}'") from None
            else:
                raise err from None
        return [v['name'] for v in tags]


    def package_download(self, pack_name: str, vers: str) -> str:
        """
        Download a package from a github repos and save it as
        <remote cache>/<organization name>/<package name>/<vers>.tar.gz

        :param str pack_name: the name of the package to download
        :param str vers: the version of the package to download
        :return: The package archive path.
        """
        url = f"{self.base_url}/repos/{self.org_name}/{pack_name}/tarball/{vers}"
        package_cache = os.path.join(self.cache, self.org_name)
        if not os.path.exists(package_cache):
            os.makedirs(package_cache)
        elif os.path.isfile(package_cache):
            raise RuntimeError(f"The tmp cache {package_cache} exists and is a file")
        tmp_archive_path = os.path.join(package_cache, f"{pack_name}-{vers}.tar.gz")
        try:
            with urllib.request.urlopen(url) as response, open(tmp_archive_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        except urllib.error.HTTPError as err:
            if 400 <= err.code < 500:
                raise RuntimeError(f"package '{pack_name}-{vers}' does not exists on repos '{self.org_name}'") \
                    from None
            else:
                raise err from None
        return tmp_archive_path


    def unarchive_package(self, path: str) -> str:
        """
        Unarchive and uncompress a package under
        <remote cache>/<organization name>/<package name>/<vers>/<package name>

        :param str path:
        :return: The path to the package
        """
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

    def __init__(self, path: str):
        """

        :param str path: The of the package root directory
        """
        self.path = os.path.realpath(path)
        self.metadata = os.path.join(self.path, 'metadata.yml')
        self.name = os.path.basename(self.path)
        self.readme = self._find_readme()
        self.check()


    def _find_readme(self) -> Any:
        """
        find the README file

        :return: The path to the README file or None if there is no file.
        """
        for ext in ('', '.md', '.rst'):
            path = os.path.join(self.path, f"README{ext}")
            if os.path.exists(path) and os.path.isfile(path):
                return path
        return None


    def _load_metadata(self) -> Dict:
        """
        Open the metadata file and de-serialize it's content
        :return:
        """
        with open(self.metadata) as raw_metadata:
            metadata = yaml.load(raw_metadata)
        return metadata


    def check(self) -> None:
        """
        Check the QA of this package
        """
        if not os.path.exists(self.path):
            raise RuntimeError()
        elif not os.path.isdir(self.path):
            raise RuntimeError()
        elif not os.path.exists(os.path.join(self.path, 'metadata.yml')):
            raise RuntimeError(f"The package '{self.name}' have no 'metadata.yml'.")
        elif not os.path.exists(os.path.join(self.path, 'definitions')):
            raise RuntimeError(f"The package '{self.name}' have no 'definitions' directory.")
        elif not os.path.isdir(os.path.join(self.path, 'definitions')):
            raise RuntimeError(f"'definitions' is not a directory.")
        elif not os.path.exists(os.path.join(self.path, 'profiles')):
            raise RuntimeError(f"The package '{self.name}' have no 'profiles' directory.")
        elif not os.path.isdir(os.path.join(self.path, 'profiles')):
            raise RuntimeError(f"'profiles' is not a directory.")
        elif not os.path.exists(os.path.join(self.path, 'LICENCE')):
            _log.warning(f"The package '{self.name}' have not any LICENCE file. May be you have not right to use it.")
        elif not self.readme:
            _log.warning(f"The package '{self.name}' have not any README file.")
        _log.info("The structure of the Package '' seems good")
        _log.info("Checking Model definitions")
        model_loc = ModelLocation(path=self.path)
        all_def = model_loc.get_all_definitions()
        model_bank = ModelBank()
        gene_bank = GeneBank()

        config = NoneConfig()
        profile_factory = ProfileFactory(config)
        parser = DefinitionParser(config, model_bank, gene_bank, profile_factory)
        parser.parse([def_loc.fqn for def_loc in all_def])
        _log.info("Definitions are consistent")


    def help(self, output=sys.stderr) -> None:
        """
        Write the content of the README file

        :param output: wher to write the help (default on stderr)
        :type output: file like object
        """
        with open(self.readme) as readme:
            for line in readme:
                print(line, file=output, end='')


    def info(self) -> str:
        """
        :return: some information about the package
        """
        metadata = self._load_metadata()
        pack_name = self.name
        cite = '\n'.join([f"\t- {c}" for c in metadata['cite']])
        info = f"""
{pack_name} {metadata['vers']}

author: {metadata['author']['name']} <{metadata['author']['email']}>

{metadata['short_desc']}

how to cite:
{cite}

documentation
\t{metadata['doc']}      

This data are realeased under {metadata['licence']}
copyright: {metadata['copyrights']}
"""
        return info


    def move(self, dest: str) -> None:
        """
        Move package from to new location *dest*
        If the destination is an existing directory,
        then the package is moved inside that directory.

        :param dest:
                    - if dest exists and is a dir, move the package to this new location
                    - if dest does not exists rename
        """
        pass