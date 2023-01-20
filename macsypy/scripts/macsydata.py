#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2023  Institut Pasteur (Paris) and CNRS.           #
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

"""
This is the entrypoint to the macsydata command
macsydata allow the user to manage the MacSyFinder models
"""

import sys
import os
import argparse
import shutil
import textwrap
import time
from typing import List, Tuple, Dict, Optional
import pathlib
import logging
import xml.etree.ElementTree as ET

import colorlog
import yaml
from packaging import requirements, specifiers, version

import macsypy
from macsypy.error import MacsyDataLimitError
from macsypy.config import MacsyDefaults, Config
from macsypy.registries import ModelRegistry, ModelLocation, scan_models_dir
from macsypy.package import RemoteModelIndex, LocalModelIndex, Package, parse_arch_path
from macsypy import licenses

# _log is set in main func
_log = None


def get_version_message():
    """
    :return: the long description of the macsyfinder version
    :rtype: str
    """
    msf_ver = macsypy.__version__
    vers_msg = f"""Macsydata {msf_ver}
Python {sys.version}

MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
See the COPYING file for details.

If you use this software please cite:
{macsypy.__citation__}
and don't forget to cite models used:
macsydata cite <model>
"""
    return vers_msg

##################
# Remote actions #
##################


def do_available(args: argparse.Namespace) -> None:
    """
    List Models available on macsy-models
    :param args: the arguments passed on the command line
    :return: None
    """
    remote = RemoteModelIndex(org=args.org)
    packages = remote.list_packages()
    for pack in packages:
        all_versions = remote.list_package_vers(pack)
        if all_versions:
            last_vers = all_versions[0]
            metadata = remote.get_metadata(pack, vers=last_vers)
            pack_vers = f"{pack} ({last_vers})"
            # 26 = length of field
            # 25 = number of displayed chars
            print(f"{pack_vers:26.25} - {metadata['short_desc']}")


def do_search(args: argparse.Namespace) -> None:
    """
    Search macsy-models for Model in a remote index.
    by default search in package name,
    if option -S is set search also in description
    by default the search is case insensitive except if
    option --match-case is set.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    try:
        remote = RemoteModelIndex(org=args.org)
        packages = remote.list_packages()
        if args.careful:
            results = _search_in_desc(args.pattern, remote, packages, match_case=args.match_case)
        else:
            results = _search_in_pack_name(args.pattern, remote, packages, match_case=args.match_case)
        for pack, last_vers, desc in results:
            pack_vers = f"{pack} ({last_vers})"
            print(f"{pack_vers:26.25} - {desc}")
    except MacsyDataLimitError as err:
        _log.critical(str(err))


def _search_in_pack_name(pattern: str, remote: RemoteModelIndex, packages: List[str],
                         match_case: bool = False) -> List[Tuple[str, str, Dict]]:
    """

    :param pattern: the substring to search packages names
    :param remote: the uri of the macsy-models index
    :param packages: list of packages to search in
    :param match_case: True if the search is case sensitive, False otherwise
    :return:
    """
    results = []
    for pack_name in packages:
        if not match_case:
            pack = pack_name.lower()
            pattern = pattern.lower()
        else:
            pack = pack_name

        if pattern in pack:
            all_versions = remote.list_package_vers(pack_name)
            if all_versions:
                metadata = remote.get_metadata(pack_name)
                last_vers = all_versions[0]
                results.append((pack_name, last_vers, metadata['short_desc']))
    return results


def _search_in_desc(pattern: str, remote: RemoteModelIndex, packages: List[str], match_case: bool = False):
    """

    :param pattern: the substring to search packages descriptions
    :param remote: the uri of the macsy-models index
    :param packages: list of packages to search in
    :param match_case: True if the search is case sensitive, False otherwise
    :return:
    """
    results = []
    for pack_name in packages:
        all_versions = remote.list_package_vers(pack_name)
        if all_versions:
            metadata = remote.get_metadata(pack_name)
            desc = metadata['short_desc']
            if not match_case:
                pack = pack_name.lower()
                desc = desc.lower()
                pattern = pattern.lower()
            else:
                pack = pack_name

            if pattern in pack or pattern in desc:
                last_vers = all_versions[0]
                results.append((pack_name, last_vers, metadata['short_desc']))
    return results


def do_download(args: argparse.Namespace) -> str:
    """
    Download tarball from remote models repository.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    try:
        remote = RemoteModelIndex(org=args.org)
        req = requirements.Requirement(args.package)
        pack_name = req.name
        specifier = req.specifier
        all_versions = remote.list_package_vers(pack_name)
        if all_versions:
            compatible_version = list(specifier.filter(all_versions))
            if compatible_version:
                vers = compatible_version[0]
                _log.info(f"Downloading {pack_name} {vers}")
                arch_path = remote.download(pack_name, vers, dest=args.dest)
                _log.info(f"Successfully downloaded packaging {pack_name} in {arch_path}")
                return arch_path
            else:
                _log.error(f"No version that satisfy requirements '{specifier}' for '{pack_name}'.")
                _log.warning(f"Available versions: {','.join(all_versions)}")
    except MacsyDataLimitError as err:
        _log.critical(str(err))


def _find_all_installed_packages(models_dir=None) -> ModelRegistry:
    """
    :return: all models installed
    """
    defaults = MacsyDefaults()
    args = argparse.Namespace()
    if models_dir is not None:
        args.models_dir = models_dir
    config = Config(defaults, args)
    model_dirs = config.models_dir()
    registry = ModelRegistry()
    for model_dir in model_dirs:
        try:
            for model_loc in scan_models_dir(model_dir, profile_suffix=config.profile_suffix()):
                registry.add(model_loc)
        except PermissionError as err:
            _log.warning(f"{model_dir} is not readable: {err} : skip it.")
    return registry


def _find_installed_package(pack_name, models_dir=None) -> Optional[ModelLocation]:
    """
    search if a package names *pack_name* is already installed

    :param pack_name: the name of the family model to search
    :return: The model location corresponding to the `pack_name`
    :rtype: :class:`macsypy.registries.ModelLocation` object
    """
    registry = _find_all_installed_packages(models_dir)
    try:
        return registry[pack_name]
    except KeyError:
        return None


def do_install(args: argparse.Namespace) -> None:
    """
    Install new models in macsyfinder local models repository.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    def clean_cache(model_index):
        if args.no_clean:
            _log.debug(f"skip cleaning {model_index.cache}")
            return
        try:
            shutil.rmtree(model_index.cache)
        except Exception:
            _log.warning(f"Cannot clean cache '{model_index.cache}': {err}")
    def create_dir(path):
        if os.path.exists(path) and not os.path.isdir(path):
            clean_cache(model_index)
            raise RuntimeError(f"'{path}' already exist and is not a directory.")
        elif not os.path.exists(path):
            os.makedirs(path)
        return path

    if os.path.exists(args.package):
        remote = False
        pack_name, inst_vers = parse_arch_path(args.package)
        user_req = requirements.Requirement(f"{pack_name}=={inst_vers}")
    else:
        remote = True
        user_req = requirements.Requirement(args.package)

    if args.target:
        dest = os.path.realpath(args.target)
        if os.path.exists(dest) and not os.path.isdir(dest):
            raise RuntimeError(f"'{dest}' already exist and is not a directory.")
        elif not os.path.exists(dest):
            os.makedirs(dest)

    pack_name = user_req.name
    inst_pack_loc = _find_installed_package(pack_name, models_dir=args.target)
    if inst_pack_loc:
        pack = Package(inst_pack_loc.path)
        try:
            local_vers = version.Version(pack.metadata['vers'])
        except FileNotFoundError:
            _log.error(f"{pack_name} locally installed is corrupted.")
            _log.warning(f"You can fix it by removing '{inst_pack_loc.path}'.")
            sys.tracebacklimit = 0
            raise RuntimeError() from None
    else:
        local_vers = None
    user_specifier = user_req.specifier
    if not user_specifier and inst_pack_loc:
        # the user do not request for a specific version
        # and there already a version installed locally
        user_specifier = specifiers.SpecifierSet(f">{local_vers}")

    if remote:
        try:
            all_available_versions = _get_remote_available_versions(pack_name, args.org)
        except (ValueError, MacsyDataLimitError) as err:
            _log.error(str(err))
            sys.tracebacklimit = 0
            raise ValueError from None
    else:
        all_available_versions = [inst_vers]

    compatible_version = list(user_specifier.filter(all_available_versions))
    if not compatible_version and local_vers:
        target_vers = version.Version(all_available_versions[0])
        if target_vers == local_vers and not args.force:
            _log.warning(f"Requirement already satisfied: {pack_name}{user_specifier} in {pack.path}.\n"
                         f"To force installation use option -f --force-reinstall.")
            return None
        elif target_vers < local_vers and not args.force:
            _log.warning(f"{pack_name} ({local_vers}) is already installed.\n"
                         f"To downgrade to {target_vers} use option -f --force-reinstall.")
            return None
        else:
            # target_vers == local_vers and args.force:
            # target_vers < local_vers and args.force:
            pass
    elif not compatible_version:
        # No compatible version and not local version
        _log.warning(f"Could not find version that satisfied '{pack_name}{user_specifier}'")
        return None
    else:
        # it exists at least one compatible version
        target_vers = version.Version(compatible_version[0])
        if inst_pack_loc:
            if target_vers > local_vers and not args.upgrade:
                _log.warning(f"{pack_name} ({local_vers}) is already installed but {target_vers} version is available.\n"
                             f"To install it please run 'macsydata install --upgrade {pack_name}'")
                return None
            elif target_vers == local_vers and not args.force:
                _log.warning(f"Requirement already satisfied: {pack_name}{user_specifier} in {pack.path}.\n"
                             f"To force installation use option -f --force-reinstall.")
                return None
            else:
                # target_vers > local_vers and args.upgrade:
                # I have to install a new package
                pass

    # if i'm here it's mean I have to install a new package
    if remote:
        _log.info(f"Downloading {pack_name} ({target_vers}).")
        model_index = RemoteModelIndex(org=args.org, cache=args.cache)
        _log.debug(f"call download with pack_name={pack_name}, vers={target_vers}")
        arch_path = model_index.download(pack_name, str(target_vers))
    else:
        model_index = LocalModelIndex(cache=args.cache)
        arch_path = args.package

    _log.info(f"Extracting {pack_name} ({target_vers}).")
    cached_pack = model_index.unarchive_package(arch_path)
    _log.debug(f"package is chached at {cached_pack}")

    if args.user:
        dest = os.path.realpath(os.path.join(os.path.expanduser('~'), '.macsyfinder', 'models'))
        create_dir(dest)
    elif args.target:
        dest = args.target
    elif 'VIRTUAL_ENV' in os.environ:
        dest = os.path.join(os.environ['VIRTUAL_ENV'], 'share', 'macsyfinder', 'models')
        create_dir(dest)
    else:
        defaults = MacsyDefaults()
        config = Config(defaults, argparse.Namespace())
        models_dirs = config.models_dir()
        if not models_dirs:
            clean_cache(model_index)
            msg = """There is no canonical directories to store models:
You can create one in your HOME to enable the models for the user 
       macsydata install --user <PACK_NAME>
or for a project 
       macsydata install --models <PACK_NAME>
In this latter case you have to specify --models-dir <path_to_models_dir> on the macsyfinder command line
for the system wide models installation please refer to the documentation.
"""
            _log.error(msg)
            sys.tracebacklimit = 0
            raise ValueError() from None
        else:
            dest = config.models_dir()[0]
    
    if inst_pack_loc:
        old_pack_path = f"{inst_pack_loc.path}.old"
        shutil.move(inst_pack_loc.path, old_pack_path)

    _log.info(f"Installing {pack_name} ({target_vers}) in {dest}")
    try:
        _log.debug(f"move {cached_pack} -> {dest}")
        shutil.move(cached_pack, dest)
    except PermissionError as err:
        clean_cache(model_index)
        _log.error(f"{dest} is not writable: {err}")
        _log.warning("Maybe you can use --user option to install in your HOME.")
        sys.tracebacklimit = 0
        raise ValueError() from None

    _log.info("Cleaning.")
    shutil.rmtree(pathlib.Path(cached_pack).parent)
    if inst_pack_loc:
        shutil.rmtree(old_pack_path)
    _log.info(f"The models {pack_name} ({target_vers}) have been installed successfully.")
    clean_cache(model_index)


def _get_remote_available_versions(pack_name, org):
    remote = RemoteModelIndex(org=org)
    all_versions = remote.list_package_vers(pack_name)
    return all_versions

#################
# Local actions #
#################


def do_uninstall(args: argparse.Namespace) -> None:
    """
    Remove models from macsyfinder local models repository.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    pack_name = args.package
    inst_pack_loc = _find_installed_package(pack_name, models_dir=args.models_dir)
    if inst_pack_loc:
        pack = Package(inst_pack_loc.path)
        shutil.rmtree(pack.path)
        _log.info(f"models '{pack_name}' in {pack.path} uninstalled.")
    else:
        _log.error(f"Models '{pack_name}' not found locally.")
        sys.tracebacklimit = 0
        raise ValueError()


def do_info(args: argparse.Namespace) -> None:
    """
    Show information about installed model.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    pack_name = args.package
    inst_pack_loc = _find_installed_package(pack_name, models_dir=args.models_dir)

    if inst_pack_loc:
        pack = Package(inst_pack_loc.path)
        print(pack.info())
    else:
        _log.error(f"Models '{pack_name}' not found locally.")
        sys.tracebacklimit = 0
        raise ValueError()


def do_list(args: argparse.Namespace) -> None:
    """
    List installed models.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    registry = _find_all_installed_packages(models_dir=args.models_dir)
    for model_loc in registry.models():
        try:
            pack = Package(model_loc.path)
            pack_vers = pack.metadata['vers']
            model_path = f"   ({model_loc.path})" if args.long else ""
            if args.outdated or args.uptodate:
                remote = RemoteModelIndex(org=args.org)
                all_versions = remote.list_package_vers(pack.name)
                specifier = specifiers.SpecifierSet(f">{pack_vers}")
                update_vers = list(specifier.filter(all_versions))
                if args.outdated and update_vers:
                    print(f"{model_loc.name}-{update_vers[0]} [{pack_vers}]{model_path}")
                if args.uptodate and not update_vers:
                    print(f"{model_loc.name}-{pack_vers}{model_path}")
            else:
                print(f"{model_loc.name}-{pack_vers}{model_path}")
        except Exception as err:
            if args.verbose > 1:
                _log.warning(str(err))


def do_freeze(args: argparse.Namespace) -> None:
    """
    display all models installed with there respective version, in requirement format.
    """
    registry = _find_all_installed_packages()
    for model_loc in sorted(registry.models(), key=lambda ml: ml.name.lower()):
        try:
            pack = Package(model_loc.path)
            pack_vers = pack.metadata['vers']
            print(f"{model_loc.name}=={pack_vers}")
        except Exception:
            pass


def do_cite(args: argparse.Namespace) -> None:
    """
    How to cite an installed model.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    pack_name = args.package
    inst_pack_loc = _find_installed_package(pack_name, models_dir=args.models_dir)
    if inst_pack_loc:
        pack = Package(inst_pack_loc.path)
        pack_citations = pack.metadata['cite']
        pack_citations = [cite.replace('\n', '\n  ') for cite in pack_citations]
        pack_citations = '\n- '.join(pack_citations)
        pack_citations = '_ ' + pack_citations.rstrip()
        macsy_cite = macsypy.__citation__
        macsy_cite = macsy_cite.replace('\n', '\n  ')
        macsy_cite = '- ' + macsy_cite
        print(f"""To cite {pack_name}:

{pack_citations}

To cite MacSyFinder:

{macsy_cite}
""")
    else:
        _log.error(f"Models '{pack_name}' not found locally.")
        sys.tracebacklimit = 0
        raise ValueError()


def do_help(args: argparse.Namespace) -> None:
    """
    Display on stdout the content of readme file
    if the readme file does nopt exists display a message to the user see :meth:`macsypy.package.help`

    :param args: the arguments passed on the command line (the package name)
    :type args: :class:`argparse.Namespace` object
    :return: None
    :raise ValueError: if the package name is not known.
    """
    pack_name = args.package
    inst_pack_loc = _find_installed_package(pack_name, models_dir=args.models_dir)
    if inst_pack_loc:
        pack = Package(inst_pack_loc.path)
        print(pack.help())
    else:
        _log.error(f"Models '{pack_name}' not found locally.")
        sys.tracebacklimit = 0
        raise ValueError()


def do_check(args: argparse.Namespace) -> None:
    """

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    pack = Package(args.path)
    errors, warnings = pack.check()
    if errors:
        for error in errors:
            _log.error(error)
        _log.error("Please fix issues above, before publishing these models.")
        sys.tracebacklimit = 0
        raise ValueError()
    if warnings:
        for warning in warnings:
            _log.warning(warning)
        _log.warning("""
macsydata says: You're only giving me a partial QA payment?
I'll take it this time, but I'm not happy.
I'll be really happy, if you fix warnings above, before to publish these models.""")

    if not warnings:
        _log.info("If everyone were like you, I'd be out of business")
        _log.info("To push the models in organization:")
        if os.path.realpath(os.getcwd()) != pack.path:
            # I use level 25 just to remove color
            _log.log(25, f"\tcd {pack.path}")
        if not os.path.exists(os.path.join(pack.path, '.git')):
            _log.info("Transform the models into a git repository")
            _log.log(25, "\tgit init .")
            _log.log(25, "\tgit add .")
            _log.log(25, "\tgit commit -m 'initial commit'")
            _log.info("add a remote repository to host the models")
            _log.info("for instance if you want to add the models to 'macsy-models'")
            _log.log(25, "\tgit remote add origin https://github.com/macsy-models/")

        _log.log(25, f"\tgit tag {pack.metadata['vers']}")
        _log.log(25, f"\tgit push origin {pack.metadata['vers']}")


def do_show_definition(args: argparse.Namespace) -> None:
    """
    display on stdout the definition if only a package or sub-package is specified
    display all model definitions in the corresponding package or subpackage

    for instance

    `TXSS+/bacterial T6SSii T6SSiii`

    display models *TXSS+/bacterial/T6SSii* and *TXSS+/bacterial/T6SSiii*

    `TXSS+/bacterial all` or `TXSS+/bacterial`

    display all models contains in *TXSS+/bacterial subpackage*

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    def display_definition(path):
        return open(path, 'r').read()

    model_family, *models = args.model
    pack_name, *sub_family = model_family.split('/')

    inst_pack_loc = _find_installed_package(pack_name, models_dir=args.models_dir)

    if inst_pack_loc:
        if not models or 'all' in models:
            root_def_name = model_family if sub_family else None
            try:
                path_2_display = sorted(
                    [(p.fqn, p.path) for p in inst_pack_loc.get_all_definitions(root_def_name=root_def_name)]
                )
            except ValueError:
                _log.error(f"'{'/'.join(sub_family)}' not found in package '{pack_name}'.")
                sys.tracebacklimit = 0
                raise ValueError() from None

            for fqn, def_path in path_2_display:
                print(f"""<!-- {fqn} {def_path} -->
{display_definition(def_path)}
""", file=sys.stdout)
        else:
            fqn_to_get = [f'{model_family}/{m}' for m in models]
            for fqn in fqn_to_get:
                try:
                    def_path = inst_pack_loc.get_definition(fqn).path
                    print(f"""<!-- {fqn} {def_path} -->
{display_definition(def_path)}
""", file=sys.stdout)
                except ValueError:
                    _log.error(f"Model '{fqn}' not found.")
                    continue
    else:
        _log.error(f"Package '{pack_name}' not found.")
        sys.tracebacklimit = 0
        raise ValueError() from None


def do_init_package(args: argparse.Namespace) -> None:
    """
    Create a template for data package

        - skeleton for metadata.yml
        - definitions directory with a skeleton of models.xml
        - profiles directory
        - skeleton for README.md file
        - COPYRIGHT file (if holders option is set)
        - LICENSE file (if license option is set)

    :param args: The parsed commandline subcommand arguments
    :return: None
    """

    def create_package_dir(package_name: str, models_dir: str = None) -> str:
        """

        :param str package_name:
        :param models_dir: the path where to create the new package
        :return: the path of the package directory
        :rtype: str
        """
        pack_path = package_name if not models_dir else os.path.join(models_dir, package_name)
        if not os.path.exists(pack_path):
            os.makedirs(pack_path)
        else:
            raise ValueError(f"{pack_path} already exist.")
        return pack_path

    def add_metadata(pack_dir: str, maintainer: str, email: str,
                     desc: str = None, license: str = None,
                     c_date: str = None, c_holders: str = None) -> None:
        """

        :param pack_dir: the package directory path
        :param maintainer: the maintainer name
        :param email: the maintainer email
        :param desc: a One line description of the package
        :param license: the license choosed
        :param c_date: the date of the copyright
        :param c_holders: the holders of the copyright
        :return: None
        """
        metadata = {
            'maintainer': {
                'name': maintainer,
                'email': email
            },
            'short_desc': desc,
            'cite': 'Place here how to cite this package, it can hold several citation',
            'doc': 'where to find documentation about this package',
            'vers': '0.1b1',
        }

        if copyright:
            metadata['copyright'] = f"Copyright (c) {c_date} {c_holders}"

        if license:
            metadata['license'] = licenses.name_2_url(license)

        with open(os.path.join(pack_dir, 'metadata.yml'), 'w') as metafile:
            yaml.dump(metadata, metafile)


    def add_def_skeleton(license: str =None) -> None:
        """
        Create a example of model definition

        :param license: the text of the license
        :return: None
        """
        model = ET.Element('model',
                           attrib={'inter_gene_max_space': "5",
                                   'min_mandatory_genes_required': "2",
                                   'min_genes_required': "3",
                                   'vers': "2.0"
                                   }
        )
        comment = ET.Comment('GENE_1 is a mandatory gene. GENE_1.hmm must exist in profiles directory')
        model.append(comment)
        mandatory = ET.SubElement(model, 'gene',
                                  attrib={'name': 'GENE_1',
                                          'presence': 'mandatory'})
        comment = ET.Comment("GENE_2 is accessory and can be exchanged with GENE_3 which play a similar role in model.\n"
                             "Both GENE_2.hmm and GENE_3.hmm must exist in profiles_directory")
        model.append(comment)
        accessory = ET.SubElement(model, 'gene',
                                  attrib={'name': 'GENE_2',
                                          'presence': 'accessory',
                                          })
        exchangeables = ET.SubElement(accessory, 'exchangeables')
        ex_gene = ET.SubElement(exchangeables, 'gene',
                                attrib={'name': 'GENE_3'})
        comment = ET.Comment("GENE_4 can be anywhere in the genome and not clusterized with some other model genes")
        model.append(comment)
        loner = ET.SubElement(model, 'gene',
                              attrib={'name': 'GENE_4',
                                      'presence': 'accessory',
                                      'loner': 'true'}
                              )
        comment = ET.Comment("GENE_5 can be shared by several systems instance from different models.")
        model.append(comment)
        multi_model = ET.SubElement(model, 'gene',
                                    attrib={'name': 'GENE_5',
                                            'presence': 'accessory',
                                            'multi_model': 'true'}
                              )
        comment = ET.Comment("GENE_6 have specific clusterisation rule")
        model.append(comment)
        inter = ET.SubElement(model, 'gene',
                              attrib={'name': 'GENE_6',
                                      'presence': 'accessory',
                                      'inter_gene_max_space': '10'}
                              )
        comment = ET.Comment("\nFor exhaustive documentation about grammar visit \n"
                             "https://macsyfinder.readthedocs.io/en/latest/modeler_guide/package.html\n")
        model.append(comment)
        tree = ET.ElementTree(model)
        try:
            ET.indent(model)
        except AttributeError:
            # workaround as ET.indent appear in python3.9
            from macsypy.utils import indent_wrapper
            ET.indent = indent_wrapper(type(tree))
            ET.indent(model)

        def_path = os.path.join(pack_dir, 'definitions', 'model_example.xml')
        tree.write(def_path,
                   encoding='UTF-8',
                   xml_declaration=True)

        if license:
            # Elementtree API does not allow to insert comment outside the tree (before root node)
            # this is the reason of this workaround
            # write the xml, read it as text, insert the comment, and write it again :-(
            with open(def_path, 'r') as def_file:
                definition = def_file.readlines()
            license = f"""<!--
{license}-->
"""
            definition.insert(1, license)
            with open(def_path, 'w') as def_path:
                def_path.writelines(definition)


    def add_license(pack_dir: str, license_text: str):
        """
        Create a license file

        :param pack_dir: the package directory path
        :param license_text: the text of the license
        :return: None
        """
        with open(os.path.join(pack_dir, 'LICENSE'), 'w') as license_file:
            license_file.write(license_text)


    def add_copyright(pack_dir: str, pack_name: str, date: str, holders: str, desc: str):
        """

        :param str pack_dir: The path of package directory
        :param str pack_name: The name of the package
        :param str date: The date (year) of package creation
        :param str holders: The copyright holders
        :param str desc: One line description of the package
        :return: None
        """
        desc = desc if desc is not None else ''
        head = textwrap.fill(f"{pack_name} - {desc}")
        text = f"""{head}
        
Copyright (c) {date} {holders}        
"""
        with open(os.path.join(pack_dir, 'COPYRIGHT'), 'w') as copyright_file:
            copyright_file.write(text)


    def add_readme(pack_dir: str, pack_name: str, desc: str):
        """

        :param str pack_dir: The path of package directory
        :param str pack_name: The name of the package
        :param str desc: One line description of the package
        :return: None
        """
        desc = desc if desc is not None else ''
        text = f"""
# {pack_name}: {desc}

Place here information about {pack_name}

- how to use it
- how to cite it
- ...

using markdown syntax
https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax
"""
        with open(os.path.join(pack_dir, 'README.md'), 'w') as readme_file:
            readme_file.write(text)


    def create_model_conf(pack_dir: str, license: str = None):
        """

        :param pack_dir: The path of the package directory
        :param license: The text of the chosen license
        :return: None
        """
        msf_defaults = MacsyDefaults()
        model_conf = ET.Element('model_config')

        weights = ET.SubElement(model_conf, 'weights')
        mandatory = ET.SubElement(weights, 'mandatory')
        mandatory.text = str(msf_defaults['mandatory_weight'])
        accessory = ET.SubElement(weights, 'accessory')
        accessory.text = str(msf_defaults['accessory_weight'])
        exchangeable = ET.SubElement(weights, 'exchangeable')
        exchangeable.text = str(msf_defaults['exchangeable_weight'])
        redundancy_penalty = ET.SubElement(weights, 'redundancy_penalty')
        redundancy_penalty.text = str(msf_defaults['redundancy_penalty'])
        out_of_cluster = ET.SubElement(weights, 'out_of_cluster')
        out_of_cluster.text = str(msf_defaults['out_of_cluster_weight'])

        filtering = ET.SubElement(model_conf, 'filtering')
        e_value_search = ET.SubElement(filtering, 'e_value_search')
        e_value_search.text = str(msf_defaults['e_value_search'])
        i_evalue_sel = ET.SubElement(filtering, 'i_evalue_sel')
        i_evalue_sel.text = str(msf_defaults['i_evalue_sel'])
        coverage_profile = ET.SubElement(filtering, 'coverage_profile')
        coverage_profile.text = str(msf_defaults['coverage_profile'])
        cut_ga = ET.SubElement(filtering, 'cut_ga')
        cut_ga.text = str(msf_defaults['cut_ga'])

        tree = ET.ElementTree(model_conf)
        conf_path = os.path.join(pack_dir, 'model_conf.xml')

        try:
            ET.indent(model_conf)
        except AttributeError:
            # workaround as ET.indent appear in python3.9
            from macsypy.utils import indent_wrapper
            ET.indent = indent_wrapper(type(tree))
            ET.indent(model_conf)

        tree.write(conf_path,
                   encoding='UTF-8',
                   xml_declaration=True)
        if license:
            # Elementtree API does not allow to insert comment outside the tree (before root node)
            # this is the reason of this workaround
            # write the xml, read it as text, insert the comment, and write it again :-(

            with open(conf_path, 'r') as conf_file:
                conf = conf_file.readlines()
            license = f"""<!--
{license}-->
"""
            conf.insert(1, license)
            with open(conf_path, 'w') as conf_file:
                conf_file.writelines(conf)


    ######################
    # Initialize Package #
    ######################
    c_date = time.localtime().tm_year
    pack_dir = create_package_dir(args.pack_name, models_dir=args.models_dir)
    def_dir = os.path.join(pack_dir, 'definitions')
    profiles_dir = os.path.join(pack_dir, 'profiles')
    license_text = None
    os.mkdir(def_dir)
    os.mkdir(profiles_dir)

    if args.holders:
        add_copyright(pack_dir, args.pack_name, c_date, args.holders, args.desc)
    else:
        _log.warning(f"Consider to add copyright to protect your rights.")

    if args.license:
        try:
            license_text = licenses.licence(args.license, args.pack_name, args.authors, c_date, args.holders, args.desc)
        except KeyError:
            _log.error(f"The license {args.license} is not managed by init (see macsydata init help). "
                       f"You will have to put the license by hand in package.")
            license_text=None
        add_license(pack_dir, license_text)
    else:
        _log.warning(f"Consider licensing {args.pack_name} to give the end-user the right to use your package,"
                     f"and protect your rights. https://data.europa.eu/elearning/en/module4/#/id/co-01")

    add_def_skeleton(license=license_text)

    create_model_conf(pack_dir, license=license_text)

    add_readme(pack_dir, args.pack_name, args.desc)

    add_metadata(pack_dir, args.maintainer, args.email, desc=args.desc, license=args.license,
                 c_date=c_date, c_holders=args.holders)

    _log.info(f"""The skeleton of {args.pack_name} is ready.
The package is located at {pack_dir}

- Edit metadata.yml and fill how to cite your package and where to find documentation about it.
- Add hmm profiles in {pack_dir}/profiles directory
- A skeleton of model definitions has been added in {pack_dir}/definitions. 
  For complete documentation about model grammar read https://macsyfinder.readthedocs.io/en/latest/modeler_guide/modeling.html
- A configuration file has been added (model_conf.xml) with default value tweak this file if needed. 
  (https://macsyfinder.readthedocs.io/en/latest/modeler_guide/package.html#model-configuration)
  
Before to publish your package you can use `macsydata check` to verify it's integrity.
"""
              )
    _log.warning("Read macsyfinder modeler guide for further details: "
                 "https://macsyfinder.readthedocs.io/en/latest/modeler_guide/index.html")

##################################
# parsing command line arguments #
##################################


def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build argument parser.

    :rtype: :class:`argparse.ArgumentParser` object
    """

    parser = argparse.ArgumentParser(
        epilog="For more details, visit the MacSyFinder website and see the MacSyFinder documentation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(r'''

         *            *               *                   * *       * 
    *           *               *   *   *  *    **                *  
      **     *    *   *  *     *                    *               *
        __  __  *         ____ *      ____    ** _  *
       |  \/  | __ _  ___/ ___| _   _|  _ | __ _| |_  __ _     *
       | |\/| |/ _` |/ __|___ \| | | | | ||/ _` |  _|/ _` |
       | |  | | (_| | (__ ___) | |_| | |_|| (_| | | | (_| |
       |_|  |_|\__,_|\___|____/ \__, |____|\__,_|_|  \__,_|
               *                |___/    *                   *
     *      *   * *     *   **         *   *  *           *
      *      *         *        *    *              *
                 *                           *  *           *     *


    MacSyData - MacSyFinder Data Management
    '''))

    # -- general options -- #

    parser.add_argument("-v", "--verbose",
                        action="count",
                        default=0,
                        help="Give more output.")
    parser.add_argument("--version",
                        action="version",
                        version=get_version_message())
    # -- subparser options -- #

    subparsers = parser.add_subparsers(help=None)
    #############
    # available #
    #############
    available_subparser = subparsers.add_parser('available',
                                                help='List Models available on macsy-models')
    available_subparser.add_argument('--org',
                                     default="macsy-models",
                                     help="The name of Model organization"
                                          "(default 'macsy-models'))"
                                     )
    available_subparser.set_defaults(func=do_available)
    ############
    # download #
    ############
    download_subparser = subparsers.add_parser('download',
                                               help='Download packages.')

    download_subparser.set_defaults(func=do_download)
    download_subparser.add_argument('-d', '--dest',
                                    default=os.getcwd(),
                                    help='Download packages into <dir>.')
    download_subparser.add_argument('--cache',
                                    help=argparse.SUPPRESS)
    download_subparser.add_argument('--org',
                                    default="macsy-models",
                                    help="The name of Model organization"
                                         "(default 'macsy-models'))"
                                    )
    download_subparser.add_argument('package', help='Package name.')
    ###########
    # Install #
    ###########
    install_subparser = subparsers.add_parser('install', help='Install packages.')
    install_subparser.set_defaults(func=do_install)
    install_subparser.add_argument('-f', '--force',
                                   action='store_true',
                                   default=False,
                                   help='Reinstall package even if it is already up-to-date.')
    install_subparser.add_argument('--org',
                                   default="macsy-models",
                                   help="The name of Model orgagnization"
                                        "(default 'macsy-models'))"
                                   )
    install_dest = install_subparser.add_mutually_exclusive_group()
    install_dest.add_argument('-u', '--user',
                              action='store_true',
                              default=False,
                              help='Install to the MacSYFinder user install directory for your platform. '
                                   'Typically ~/.macsyfinder/data')
    install_dest.add_argument('-t', '--target', '--models-dir',
                              dest='target',
                              help='Install packages into <TARGET> dir instead in canonical location')

    install_subparser.add_argument('-U', '--upgrade',
                                   action='store_true',
                                   default=False,
                                   help='Upgrade specified package to the newest available version.')
    install_subparser.add_argument('package',
                                   help='Package name.')
    install_subparser.add_argument('--cache',
                                   help=argparse.SUPPRESS)
    install_subparser.add_argument('--no-clean',
                                   action='store_true',
                                   default=False,
                                   # do not clean cache for debugging purpose ONLY
                                   help=argparse.SUPPRESS)
    #############
    # Uninstall #
    #############
    uninstall_subparser = subparsers.add_parser('uninstall',
                                                help='Uninstall packages.')
    uninstall_subparser.set_defaults(func=do_uninstall)
    uninstall_subparser.add_argument('package',
                                     help='Package name.')
    uninstall_subparser.add_argument('--target, --models-dir',
                                     dest='models_dir',
                                     help='the path of the alternative root directory containing package instead used '
                                     'canonical locations')
    ##########
    # search #
    ##########
    search_subparser = subparsers.add_parser('search',
                                             help='Discover new packages.')
    search_subparser.set_defaults(func=do_search)
    search_subparser.add_argument('--org',
                                  default="macsy-models",
                                  help="The name of Model organization"
                                       "(default macsy-models))"
                                  )
    search_subparser.add_argument('-S', '--careful',
                                  default=False,
                                  action='store_true',
                                  help='')
    search_subparser.add_argument('--match-case',
                                  default=False,
                                  action='store_true',
                                  help='')
    search_subparser.add_argument('pattern',
                                  help='Searches for packages matching the pattern.')
    ########
    # info #
    ########
    info_subparser = subparsers.add_parser('info',
                                           help='Show information about packages.')
    info_subparser.add_argument('package',
                                help='Package name.')
    info_subparser.set_defaults(func=do_info)
    info_subparser.add_argument('--models-dir',
                                help='the path of the alternative root directory containing package instead used '
                                     'canonical locations')
    ########
    # list #
    ########
    list_subparser = subparsers.add_parser('list',
                                           help='List installed packages.')
    list_subparser.set_defaults(func=do_list)
    list_subparser.add_argument('-o', '--outdated',
                                action='store_true',
                                default=False,
                                help='List outdated packages.')
    list_subparser.add_argument('-u', '--uptodate',
                                action='store_true',
                                default=False,
                                help='List uptodate packages')
    list_subparser.add_argument('--org',
                                default="macsy-models",
                                help="The name of Model organization"
                                     "(default macsy-models))"
                                )
    list_subparser.add_argument('--models-dir',
                                help='the path of the alternative root directory containing package instead used '
                                     'canonical locations')
    list_subparser.add_argument('--long', '-l',
                                action='store_true',
                                default=False,
                                help="in addition displays the path where is store each package"
                                )
    list_subparser.add_argument('-v',
                                dest='long',
                                action='store_true',
                                default=False,
                                help="alias for -l/--long option"
                                )
    ##########
    # freeze #
    ##########
    freeze_subparser = subparsers.add_parser('freeze',
                                             help='List installed models in requirements format.')
    freeze_subparser.add_argument('--models-dir',
                                   help='the path of the alternative root directory containing package instead used '
                                        'canonical locations')
    freeze_subparser.set_defaults(func=do_freeze)
    ########
    # cite #
    ########
    cite_subparser = subparsers.add_parser('cite',
                                           help='How to cite a package.')
    cite_subparser.set_defaults(func=do_cite)
    cite_subparser.add_argument('--models-dir',
                                help='the path of the alternative root directory containing package instead used '
                                     'canonical locations')
    cite_subparser.add_argument('package',
                                help='Package name.')
    ########
    # help #
    ########
    help_subparser = subparsers.add_parser('help',
                                           help='get online documentation.')
    help_subparser.set_defaults(func=do_help)
    help_subparser.add_argument('package',
                                help='Package name.')
    help_subparser.add_argument('--models-dir',
                                help='the path of the alternative root directory containing package instead used '
                                     'canonical locations')
    #########
    # check #
    #########
    check_subparser = subparsers.add_parser('check',
                                            help='check if the directory is ready to be publish as data package')
    check_subparser.set_defaults(func=do_check)
    check_subparser.add_argument('path',
                                 nargs='?',
                                 default=os.getcwd(),
                                 help='the path to root directory models to check')

    ##############
    # definition #
    ##############
    def_subparser = subparsers.add_parser('definition',
                                            help='show a model definition ')
    def_subparser.set_defaults(func=do_show_definition)
    def_subparser.add_argument('model',
                               nargs='+',
                               help='the family and name(s) of a model(s) eg: TXSS T6SS T4SS or TFF/bacterial T2SS')
    def_subparser.add_argument('--models-dir',
                               help='the path to the alternative root directory containing packages instead to the '
                                    'canonical locations')
    ########
    # init #
    ########
    init_subparser = subparsers.add_parser('init',
                                           help='Create a template for a new data package')
    init_subparser.set_defaults(func=do_init_package)
    init_subparser.add_argument('--pack-name',
                                required=True,
                                help='The name of the data package.')
    init_subparser.add_argument('--maintainer',
                                required=True,
                                help='The name of the package maintainer.')
    init_subparser.add_argument('--email',
                                required=True,
                                help='The email of the package maintainer.')
    init_subparser.add_argument('--authors',
                                required=True,
                                help="The authors of the package. Could be different that the maintainer."
                                     "Could be several persons. Surround the names by quotes 'John Doe, Richard Miles'")
    init_subparser.add_argument('--license',
                                choices=['cc-by', 'cc-by-sa', 'cc-by-nc', 'cc-by-nc-sa', 'cc-by-nc-nd'],
                                help="""The license under this work will be released.
if the license you choice is not in the list, you can do it manually
by adding the license file in package and add suitable headers in model definitions.""")
    init_subparser.add_argument('--holders',
                                help="The holders of the copyright")
    init_subparser.add_argument('--desc',
                                help="A short description (one line) of the package")
    init_subparser.add_argument('--models-dir',
                                help='The path of an alternative models directory by default the package will be created here.' )
    return parser


def cmd_name(args: argparse.Namespace) -> str:
    """
    Return the name of the command being executed
    (scriptname + operation).

    Example
        macsydata uninstall

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: str
    """
    assert 'func' in args
    func_name = args.func.__name__.replace('do_', '')
    return f"macsydata {func_name}"


def init_logger(level='INFO', out=True):
    """

    :param level: The logger threshold could be a positive int or string
                  among: 'CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'
    :param out: if the log message must be displayed
    :return: logger
    :rtype: :class:`logging.Logger` instance
    """

    logger = colorlog.getLogger('macsydata')
    handlers = []
    if out:
        stdout_handler = colorlog.StreamHandler(sys.stderr)
        if level <= logging.DEBUG:
            msg_formatter = "%(log_color)s%(levelname)-8s : %(module)s: L %(lineno)d :%(reset)s %(message)s"
        else:
            msg_formatter = "%(log_color)s%(message)s"
        stdout_formatter = colorlog.ColoredFormatter(msg_formatter,
                                                     datefmt=None,
                                                     reset=True,
                                                     log_colors={
                                                         'DEBUG': 'cyan',
                                                         'INFO': 'green',
                                                         'WARNING': 'yellow',
                                                         'ERROR': 'red',
                                                         'CRITICAL': 'bold_red',
                                                     },
                                                     secondary_log_colors={},
                                                     style='%'
                                                     )
        stdout_handler.setFormatter(stdout_formatter)
        logger.addHandler(stdout_handler)
        handlers.append(stdout_handler)
    else:
        null_handler = logging.NullHandler()
        logger.addHandler(null_handler)
        handlers.append(null_handler)
    if isinstance(level, str):
        level = getattr(logging, level)
    logger.setLevel(level)
    return logger


def verbosity_to_log_level(verbosity: int) -> int:
    """
    transform the number of -v option in loglevel
    :param int verbosity: number of -v option on the command line
    :return: an int corresponding to a logging level
    """
    level = max((logging.INFO - (10 * verbosity), 1))
    return level


def main(args=None) -> None:
    """
    Main entry point.

    :param args: the arguments passed on the command line (before parsing)
    :type args: list
    :rtype: int
    """
    global _log
    args = sys.argv[1:] if args is None else args
    parser = build_arg_parser()
    parsed_args = parser.parse_args(args)
    log_level = verbosity_to_log_level(parsed_args.verbose)
    # set logger for module 'package'
    macsypy.init_logger()
    macsypy.logger_set_level(level=log_level)
    # set logger for this script
    _log = init_logger(log_level)

    if 'func' in parsed_args:
        parsed_args.func(parsed_args)
        _log.debug(f"'{cmd_name(parsed_args)}' command completed successfully.")
    else:
        # macsydata command is run without any subcommand
        parser.print_help()


if __name__ == "__main__":
    main()
