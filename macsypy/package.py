#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2021  Institut Pasteur (Paris) and CNRS.           #
# See the COPYRIGHT file for details                                    #
#                                                                       #
# This file is part of MacSyFinder package.                             #
#                                                                       #
# MacSyFinder is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# MacSyFinder is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details .                         #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with MacSyFinder (COPYING).                                     #
# If not, see <https://www.gnu.org/licenses/>.                          #
#########################################################################

import os
import tempfile
import urllib.request
import urllib.parse
import json
import yaml
import shutil
import tarfile
import copy
import abc
from typing import List, Dict, Tuple, Optional
import logging
_log = logging.getLogger(__name__)

from .config import NoneConfig
from .registries import ModelLocation, ModelRegistry
from .profile import ProfileFactory
from .definition_parser import DefinitionParser
from .model import ModelBank
from .gene import GeneBank
from .model_conf_parser import ModelConfParser
from .error import MacsydataError, MacsyDataLimitError, MacsypyError


class AbstractModelIndex(metaclass=abc.ABCMeta):
    """
    This the base class for ModelIndex.
    This class cannot be implemented, it must be subclassed
    """

    def __new__(cls, *args, **kwargs):
        if cls.__bases__ == (object,):
            raise TypeError(f'{cls.__name__} is abstract cannot be instantiated.')
        return super(AbstractModelIndex, cls).__new__(cls)


    def __init__(self, cache: str = ''):
        """

        """
        self.org_name: str = None
        if cache:
            self.cache: str = cache
        else:
            self.cache = os.path.join(tempfile.gettempdir(), 'tmp-macsy-cache')


    def unarchive_package(self, path: str) -> str:
        """
        Unarchive and uncompress a package under
        `<remote cache>/<organization name>/<package name>/<vers>/<package name>`

        :param str path:
        :return: The path to the package
        """
        name, vers = parse_arch_path(path)
        dest_dir = os.path.join(self.cache, self.org_name, name, vers)
        dest_unarchive_path = os.path.join(dest_dir, name)
        if os.path.exists(dest_unarchive_path):
            _log.info(f"Removing old models {dest_unarchive_path}")
            shutil.rmtree(dest_unarchive_path)
        tar = tarfile.open(path, 'r:gz')
        tar_dir_name = tar.next().name
        tar.extractall(path=dest_dir)
        # github prefix the archive root directory with the organization name
        # add suffix with a random suffix
        # for instance for TXSS models
        # the unarchive will named macsy-models-TXSS-64889bd
        unarchive_pack = os.path.join(dest_dir, tar_dir_name)
        if unarchive_pack != dest_unarchive_path:
            os.rename(unarchive_pack, dest_unarchive_path)
        return dest_unarchive_path


class LocalModelIndex(AbstractModelIndex):
    """
    It allow to manage installation from a local package (tarball)
    """

    def __init__(self, cache=None) -> None:
        """

        """
        super().__init__(cache=cache)
        self.org_name: str = 'local'


class RemoteModelIndex(AbstractModelIndex):
    """
    This class allow to interact with ModelIndex on github
    """

    def __init__(self, org: str = "macsy-models", cache=None) -> None:
        """

        :param org: The name of the organization on github where are stored the models
        """
        super().__init__(cache=cache)
        self.org_name = urllib.parse.quote(org)
        self.base_url: str = "https://api.github.com"
        if not self.remote_exists():
            raise ValueError(f"the '{self.org_name}' organization does not exist.")


    def _url_json(self, url: str) -> Dict:
        """
        Get the url, deserialize the data as json

        :param str url: the url to download
        :return: the json corresponding to the response url
        """
        try:
            r = urllib.request.urlopen(url).read()
        except urllib.error.HTTPError as err:
            if err.code == 403:
                raise MacsyDataLimitError("You reach the maximum number of request per hour to github.\n"
                                          "Please wait before to try again.") from None
            else:
                raise err
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
            vers = versions[0]
        else:
            if vers not in versions:
                raise RuntimeError(f"The version '{vers}' does not exists for model {pack_name}.")
        pack_name = urllib.parse.quote(pack_name)
        vers = urllib.parse.quote(vers)
        metadata_url = f"https://raw.githubusercontent.com/{self.org_name}/{pack_name}/{vers}/metadata.yml"
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
        pack_name = urllib.parse.quote(pack_name)
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
        _log.debug(f"call download with pack_name={pack_name}, vers={vers}, dest={dest}")
        safe_pack_name = urllib.parse.quote(pack_name)
        safe_vers = urllib.parse.quote(vers)
        url = f"{self.base_url}/repos/{self.org_name}/{safe_pack_name}/tarball/{safe_vers}"
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


class Package:
    """
    This class Modelize a package of Models
    a package is a directory with the name of the models family
    it must contains at least
    - a subdirectory definitions
    - a subdirectory profiles
    - a file metadata.yml
    it is also recomanded to add a file
    for licensing and copyright and a README.
    for further explanation see TODO

    """

    def __init__(self, path: str) -> None:
        """

        :param str path: The of the package root directory
        """
        self.path: str = os.path.realpath(path)
        self.metadata_path: str = os.path.join(self.path, 'metadata.yml')
        self._metadata: Dict = None
        self.name: str = os.path.basename(self.path)
        self.readme: str = self._find_readme()


    def _find_readme(self) -> Optional[str]:
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
        """

        :return: The parsed metadata as a dict
        """
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


    def check(self) -> Tuple[List[str], List[str]]:
        """
        Check the QA of this package
        """
        all_warnings = []
        all_errors = []
        for meth in self._check_structure, self._check_metadata, self._check_model_consistency, self._check_model_conf:
            errors, warnings = meth()
            all_errors.extend(errors)
            all_warnings.extend(warnings)
            if all_errors:
                break
        return all_errors, all_warnings


    def _check_structure(self) -> Tuple[List[str], List[str]]:
        """
        Check the QA structure of the package

        :return: errors and warnings
        :rtype: tuple of 2 lists ([str error_1, ...], [str warning_1, ...])
        """
        _log.info(f"Checking '{self.name}' package structure")
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

            if not os.path.exists(os.path.join(self.path, 'LICENSE')):
                warnings.append(f"The package '{self.name}' have not any LICENSE file. "
                                f"May be you have not right to use it.")
            if not self.readme:
                warnings.append(f"The package '{self.name}' have not any README file.")
        return errors, warnings


    def _check_model_consistency(self) -> Tuple[List, List]:
        """
        check if each xml seems well write, each genes have an associated profile, etc

        :return:
        """
        _log.info(f"Checking '{self.name}' Model definitions")
        model_loc = ModelLocation(path=self.path)
        all_def = model_loc.get_all_definitions()
        model_bank = ModelBank()
        gene_bank = GeneBank()

        config = NoneConfig()
        config.models_dir = lambda: self.path
        try:
            profile_factory = ProfileFactory(config)
            model_registry = ModelRegistry()
            model_registry.add(model_loc)
            parser = DefinitionParser(config, model_bank, gene_bank, model_registry, profile_factory)
            parser.parse(all_def)
        finally:
            del config.models_dir
        _log.info("Definitions are consistent")
        # to respect same api as _check_metadata and _check_structure
        return [], []


    def _check_model_conf(self) -> Tuple[List[str], List[str]]:
        """
        check if a model configuration file is present in the package (model_conf.xml)
        if the syntax of this file is good.

        :return:
        """
        _log.info(f"Checking '{self.name}' model configuration")
        errors = []
        warnings = []
        conf_file = os.path.join(self.path, 'model_conf.xml')
        if os.path.exists(conf_file):
            mcp = ModelConfParser(conf_file)
            try:
                mcp.parse()
            except (ValueError, MacsypyError) as err:
                errors.append(str(err))
        else:
            _log.info(f"There is no model configuration for package {self.name}.")
        return errors, warnings


    def _check_metadata(self) -> Tuple[List[str], List[str]]:
        """
        Check the QA of package metadata_path

        :return: errors and warnings
        :rtype: tuple of 2 lists ([str error_1, ...], [str warning_1, ...])
        """
        _log.info(f"Checking '{self.name}' {self.metadata_path}")
        errors = []
        warnings = []
        data = self._load_metadata()
        must_have = ("maintainer", "short_desc", "vers")
        nice_to_have = ("cite", "doc", "license", "copyright")
        for item in must_have:
            if item not in data:
                errors.append(f"field '{item}' is mandatory in {self.metadata_path}.")
        for item in nice_to_have:
            if item not in data:
                warnings.append(f"It's better if the field '{item}' is setup in {self.metadata_path} file")
        if "maintainer" in data:
            for item in ("name", "email"):
                if item not in data["maintainer"]:
                    errors.append(f"field 'maintainer.{item}' is mandatory in {self.metadata_path}.")
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
        if 'license' not in metadata:
            metadata['license'] = "No license available"
        copyrights = f"copyright: {metadata['copyright']}" if 'copyright' in metadata else ''
        pack_name = self.name
        cite = '\n'.join([f"\t- {c}".replace('\n', '\n\t  ') for c in metadata['cite']]).rstrip()
        info = f"""
{pack_name} ({metadata['vers']})

maintainer: {metadata['maintainer']['name']} <{metadata['maintainer']['email']}>

{metadata['short_desc']}

how to cite:
{cite}

documentation
\t{metadata['doc']}

This data are released under {metadata['license']}
{copyrights}
"""
        return info


def parse_arch_path(path: str) -> Tuple[str, str]:
    """

    :param str path: the path to the archive
    :return: the name of the package and it's version
    :rtype: tuple
    :raise ValueError: if the extension of the package is neither '.tar.gz' nor '.tgz'
                       or if the package does not seem to include version 'pack_name-<vers>.ext'
    """
    pack_vers_name = os.path.basename(path)
    if pack_vers_name.endswith('.tar.gz'):
        pack_vers_name = pack_vers_name[:-7]
    elif pack_vers_name.endswith('.tgz'):
        pack_vers_name = pack_vers_name[:-4]
    else:
        raise ValueError(f"{path} does not seem to be a package (a tarball).")
    *pack_name, vers = pack_vers_name.split('-')
    if not pack_name:
        raise ValueError(f"{path} does not seem to not be versioned.")
    pack_name = '-'.join(pack_name)
    return pack_name, vers
