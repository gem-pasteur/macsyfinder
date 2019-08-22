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
import copy
from typing import List, Dict, Any

import logging
_log = logging.getLogger(__name__)

from .config import NoneConfig
from .registries import ModelLocation, ModelRegistry
from .definition_parser import DefinitionParser
from .model import ModelBank
from .gene import GeneBank, ProfileFactory
from .error import MacsydataError


class RemoteModelIndex:

    def __init__(self, org: str = "macsy-models"):
        """

        :param org: The name of the organization on github where are stored the models
        """
        self.org_name = org
        self.base_url = "https://api.github.com"
        self.cache = os.path.join(tempfile.gettempdir(), 'tmp-macsy-cache')
        if not self.remote_exists():
            raise ValueError(f"the '{self.org_name}' organization does not exist.")


    def _url_json(self, url: str) -> Dict:
        """
        Get the url, deserialize the data as json

        :param str url: the url to dowload
        :return: the json corresponding to the response url
        """
        r = urllib.request.urlopen(url).read()
        j = json.loads(r.decode('utf-8'))
        return j


    def remote_exists(self) -> bool:
        """
        check if the remote exists and is an organization
        :return: True if the Remote url point to a github Organization, False otherwise
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

    def get_metadata(self, pack_name: str, vers: str = 'latest') -> Dict:
        """
        Fetch the metadata_path from a remote package

        :param str pack_name: The package name
        :param str vers: The package version
        :return: the metadata_path corresponding to this package/version
        :rtype: dictionary corresponding of the yaml parsing of the metadata_path file.
        """
        versions = self.list_package_vers(pack_name)
        if not versions:
            raise MacsydataError(f"No official version available for model '{pack_name}'")
        elif vers == 'latest':
            vers == versions[0]
        else:
            if vers not in versions:
                raise RuntimeError(f"The version '{vers}' does not exists for model {pack_name}.")
        metadata_url = f"https://raw.githubusercontent.com/{self.org_name}/{pack_name}/{versions[0]}/metadata.yml"
        try:
            with urllib.request.urlopen(metadata_url) as response:
                metadata = response.read().decode("utf-8")
        except urllib.error.HTTPError as err:
            if 400 < err.code < 500:
                raise MacsydataError(f"cannot fetch '{metadata_url}' check '{pack_name}'")
            elif err.code >= 500:
                raise err from None
            else:
                raise err from None
        metadata = yaml.safe_load(metadata)
        return metadata


    def list_packages(self) -> List[str]:
        """
        list all model packages available on a model repos
        :return: The list of package names.
        """
        url = f"{self.base_url}/orgs/{self.org_name}/repos"
        _log.debug(f"get {url}")
        packages = self._url_json(url)
        return [p['name'] for p in packages]


    def list_package_vers(self, pack_name: str) -> List[str]:
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
                raise ValueError(f"package '{pack_name}' does not exists on repos '{self.org_name}'") from None
            else:
                raise err from None
        return [v['name'] for v in tags]


    def download(self, pack_name: str, vers: str, dest: str = None) -> str:
        """
        Download a package from a github repos and save it as
        <remote cache>/<organization name>/<package name>/<vers>.tar.gz

        :param str pack_name: the name of the package to download
        :param str vers: the version of the package to download
        :param str dest: The path to the directory where save the package
                         This directory must exists
                         If dest is None, the macsyfinder cache will be used
        :return: The package archive path.
        """
        url = f"{self.base_url}/repos/{self.org_name}/{pack_name}/tarball/{vers}"
        if not dest:
            package_cache = os.path.join(self.cache, self.org_name)
            if os.path.exists(self.cache) and not os.path.isdir(self.cache):
                raise NotADirectoryError(f"The tmp cache '{self.cache}' already exists")
            elif not os.path.exists(package_cache):
                os.makedirs(package_cache)
            tmp_archive_path = os.path.join(package_cache, f"{pack_name}-{vers}.tar.gz")
        else:
            tmp_archive_path = os.path.join(dest, f"{pack_name}-{vers}.tar.gz")
        try:
            with urllib.request.urlopen(url) as response, open(tmp_archive_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        except urllib.error.HTTPError as err:
            if 400 <= err.code < 500:
                raise ValueError(f"package '{pack_name}-{vers}' does not exists on repos '{self.org_name}'") \
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
        *name, vers = '.'.join(os.path.basename(path).split('.')[:-2]).split('-')
        name = '-'.join(name)
        dest_dir = os.path.join(self.cache, self.org_name, name, vers)
        tar = tarfile.open(path, 'r|gz')
        tar.extractall(path=dest_dir)
        src = glob.glob(os.path.join(dest_dir, f"{self.org_name}-{name}-*"))
        if len(src) == 1:
            src = src[0]
        elif len(src) > 1:
            raise MacsydataError(f"Too many matching packages. May be you have to clean {dest_dir}")
        else:
            raise MacsydataError("An error occurred during archive extraction")
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
        self.metadata_path = os.path.join(self.path, 'metadata.yml')
        self._metadata = None
        self.name = os.path.basename(self.path)
        self.readme = self._find_readme()


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

    @property
    def metadata(self) -> Dict:
        if not self._metadata:
            self._metadata = self._load_metadata()
        # to avoid side effect
        return copy.deepcopy(self._metadata)


    def _load_metadata(self) -> Dict:
        """
        Open the metadata_path file and de-serialize it's content
        :return:
        """
        with open(self.metadata_path) as raw_metadata:
            metadata = yaml.safe_load(raw_metadata)
        return metadata


    def check(self) -> List:
        """
        Check the QA of this package
        """
        all_warnings = []
        all_errors = []
        for meth in self._check_structure, self._check_metadata, self._check_model_consistency:
            errors, warnings = meth()
            all_errors.extend(errors)
            all_warnings.extend(warnings)
            if all_errors:
                break
        return all_errors, all_warnings


    def _check_structure(self) -> List[str]:
        """
        Check the QA structure of the package

        :return: errors and warnings
        :rtype: tuple of 2 lists ([str error_1, ...], [str warning_1, ...])
        """
        _log.info(f"Checking '{self.name}'package structure")
        errors = []
        warnings = []
        if not os.path.exists(self.path):
            errors.append(f"The package '{self.name}' does not exists.")
        elif not os.path.isdir(self.path):
            errors.append(f"'{self.name}' is not a directory ")
        elif not os.path.exists(os.path.join(self.path, 'metadata.yml')):
            errors.append(f"The package '{self.name}' have no 'metadata.yml'.")
        if not errors:
            # check several criteria and don't stop at the first problem.
            # this is why I use several If and not one set of if/elif
            if not os.path.exists(os.path.join(self.path, 'definitions')):
                errors.append(f"The package '{self.name}' have no 'definitions' directory.")
            elif not os.path.isdir(os.path.join(self.path, 'definitions')):
                errors.append(f"'{os.path.join(self.path, 'definitions')}' is not a directory.")

            if not os.path.exists(os.path.join(self.path, 'profiles')):
                errors.append(f"The package '{self.name}' have no 'profiles' directory.")
            elif not os.path.isdir(os.path.join(self.path, 'profiles')):
                errors.append(f"'{os.path.join(self.path, 'profiles')}' is not a directory.")

            if not os.path.exists(os.path.join(self.path, 'LICENCE')):
                warnings.append(f"The package '{self.name}' have not any LICENCE file. "
                                f"May be you have not right to use it.")
            if not self.readme:
                warnings.append(f"The package '{self.name}' have not any README file.")
        return errors, warnings


    def _check_model_consistency(self) -> None:
        _log.info(f"Checking '{self.name}' Model definitions")
        model_loc = ModelLocation(path=self.path)
        all_def = model_loc.get_all_definitions()
        model_bank = ModelBank()
        gene_bank = GeneBank()

        config = NoneConfig()
        config.models_dir = lambda: self.path
        profile_factory = ProfileFactory(config)
        model_registry = ModelRegistry()
        model_registry.add(model_loc)
        parser = DefinitionParser(config, model_bank, gene_bank, profile_factory, model_registry)
        parser.parse([def_loc.fqn for def_loc in all_def])
        _log.info("Definitions are consistent")
        # to respect same api as _check_metadata and _check_structure
        return [], []


    def _check_metadata(self) -> List[str]:
        """
        Check the QA of package metadata_path

        :return: errors and warnings
        :rtype: tuple of 2 lists ([str error_1, ...], [str warning_1, ...])
        """
        _log.info(f"Checking '{self.name}' metadata_path")
        errors = []
        warnings = []
        data = self._load_metadata()
        must_have = ("author", "short_desc", "vers" )
        nice_to_have = ("cite", "doc", "licence", "copyright")
        for item in must_have:
            if item not in data:
                errors.append(f"field '{item}' is mandatory in metadata_path.")
        for item in nice_to_have:
            if item not in data:
                warnings.append(f"It's better if the field '{item}' is setup in metadata_path file")
        if "author" in data:
            for item in ("name", "email"):
                if item not in data["author"]:
                    errors.append(f"field 'author.{item}' is mandatory in metadata_path.")
        return errors, warnings


    def help(self) -> str:
        """
        return the content of the README file
        """
        if self.readme:
            with open(self.readme) as readme:
                pack_help = ''.join(readme.readlines())
        else:
            pack_help = f"No help available for package '{self.name}'."
        return pack_help


    def info(self) -> str:
        """
        :return: some information about the package
        """
        metadata = self._load_metadata()
        if 'cite' not in metadata:
            metadata['cite'] = ["No citation available\n"]
        if 'doc' not in metadata:
            metadata['doc'] = "No documentation available"
        if 'licence' not in metadata:
            metadata['licence'] = "No licence available"
        copyrights = f"copyright: {metadata['copyright']}" if 'copyright' in metadata else ''
        pack_name = self.name
        cite = '\n'.join([f"\t- {c}".replace('\n', '\n\t  ') for c in metadata['cite']])
        info = f"""
{pack_name} ({metadata['vers']})

author: {metadata['author']['name']} <{metadata['author']['email']}>

{metadata['short_desc']}

how to cite:
{cite}
documentation
\t{metadata['doc']}

This data are released under {metadata['licence']}
{copyrights}
"""
        return info
