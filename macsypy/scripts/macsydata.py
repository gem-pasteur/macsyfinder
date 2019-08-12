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


import sys
import os
import argparse
import logging
import shutil
from textwrap import dedent
from typing import List
import pathlib

import colorlog
_log = colorlog.getLogger('macsypy')
from packaging import requirements
from packaging import specifiers

import macsypy
from macsypy.config import MacsyDefaults, Config
from macsypy.registries import ModelRegistry , scan_models_dir
from macsypy.package import RemoteModelIndex, Package


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
    remote = RemoteModelIndex(org=args.org)
    packages = remote.list_packages()
    if args.careful:
        results = _search_in_desc(args.pattern, remote, packages, match_case=args.match_case)
    else:
        results = _search_in_pack_name(args.pattern, remote, packages, match_case=args.match_case)
    for pack, last_vers, desc in results:
        pack_vers = f"{pack} ({last_vers})"
        print(f"{pack_vers:26.25} - {desc}")


def _search_in_pack_name(pattern: str, remote: RemoteModelIndex, packages: List[str], match_case: bool = False) -> List[str]:
    """

    :param pattern:
    :param remote:
    :param packages:
    :param match_case:
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
    remote = RemoteModelIndex(org=args.org)
    req = requirements.Requirement(args.package)
    pack_name = req.name
    specifier = req.specifier
    all_versions = remote.list_package_vers(pack_name)
    if all_versions:
        compatible_version = list(specifier.filter(all_versions))
        if compatible_version:
            vers = compatible_version[0]
            print(f"Downloading {pack_name} {vers}")
            arch_path = remote.download(pack_name, vers, dest=args.dest)
            print(f"Successfully downloaded packaging {pack_name} in {arch_path}")
            return arch_path
        else:
            raise ValueError(f"No version that satisfy requirements '{specifier}' for '{pack_name}'. "
                             f"Available versions: {','.join(all_versions)}")


def _find_local_package(pack_name):
    defaults = MacsyDefaults()
    config = Config(defaults, argparse.Namespace())
    system_model_dir = config.models_dir()
    user_model_dir = os.path.join(os.path.expanduser('~'), '.macsyfinder', 'data')
    model_dirs = (system_model_dir, user_model_dir) if os.path.exists(user_model_dir) else (system_model_dir,)
    registry = ModelRegistry()
    for model_dir in model_dirs:
        for model_loc in scan_models_dir(model_dir, profile_suffix=config.profile_suffix):
            registry.add(model_loc)
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
    req = requirements.Requirement(args.package)
    pack_name = req.name
    specifier = req.specifier
    local_pack = _find_local_package(pack_name)

    if local_pack:
        pack = Package(local_pack.path, check=False)
        local_vers = pack.metadata['vers']

    remote = RemoteModelIndex(org=args.org)
    all_versions = remote.list_package_vers(pack_name)
    if not specifier and local_pack:
        specifier = specifiers.SpecifierSet(f">{local_vers}")
    compatible_version = list(specifier.filter(all_versions))
    if local_pack and compatible_version and not args.upgrade:
        print(f"{pack_name} is already installed but {compatible_version[0]} is available.\n"
              f"To install it please run 'macsydata install --upgrade {pack_name}'")
        return None
    elif local_pack and not compatible_version:
        if args.force:
            print("#### ", args.force)
            if local_vers in all_versions:
                target_vers = local_vers
            elif all_versions:
                target_vers = all_versions[0]
            else:
                print(f"Cannot find version for models '{pack_name}'.")
                return None
        else:
            print(f"Requirement already satisfied: {pack_name}{specifier} in {pack.path}.\n"
                  f"To force installation use option -f --force-reinstall.")
            return None
    elif not local_pack and not compatible_version:
        raise RuntimeError(f"Could not find version that satisfied '{pack_name}{specifier}'")
    else:
        target_vers = compatible_version[0]

    print(f"Downloading {pack_name} ({target_vers}).")
    tmp_arch = remote.download(pack_name, target_vers)
    print(f"Extracting {pack_name} ({target_vers}).")
    cached_pack = remote.unarchive_package(tmp_arch)
    if args.user:
        dest = os.path.join(os.path.expanduser('~'), '.macsyfinder', 'data')
    else:
        defaults = MacsyDefaults()
        config = Config(defaults, argparse.Namespace())
        dest = config.models_dir()
    if local_pack:
        old_pack_path = f"{local_pack.path}.old"
        shutil.move(local_pack.path, old_pack_path)
    print(f"Installing {pack_name} ({target_vers}) in {dest}")
    shutil.move(cached_pack, dest)
    if local_pack:
        print("Cleaning.")
        shutil.rmtree(old_pack_path)
        shutil.rmtree(pathlib.Path(cached_pack).parent)
    print(f"The models {pack_name} ({target_vers}) have been installed successfully.")


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
    local_pack = _find_local_package(pack_name)

    if local_pack:
        pack = Package(local_pack.path, check=False)
        shutil.rmtree(pack.path)
        print(f"models '{pack_name}' in {pack.path} uninstalled.")
    else:
        print(f"Models '{pack_name}' not found locally.")


def do_info(args: argparse.Namespace) -> None:
    """
    Show information about installed model.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    pack_name = args.package
    local_pack = _find_local_package(pack_name)

    if local_pack:
        pack = Package(local_pack.path, check=False)
        print(pack.info())
    else:
        print(f"Models '{pack_name}' not found locally.")


def do_list(args: argparse.Namespace) -> None:
    """
    List installed models.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    raise Exception('Not implemented')


def do_cite(args: argparse.Namespace) -> None:
    """
    How to cite an installed model.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    pack_name = args.package
    local_pack = _find_local_package(pack_name)
    if local_pack:
        pack = Package(local_pack.path, check=False)
        pack_citations = pack.metadata['cite']
        pack_citations = [cite.replace('\n', '\n  ') for cite in pack_citations]
        pack_citations = '\n- '.join(pack_citations)
        pack_citations = '_ ' + pack_citations
        macsy_cite = macsypy.__citation__
        macsy_cite = macsy_cite.replace('\n', '\n  ')
        macsy_cite = '- ' + macsy_cite
        print(f"""To cite {pack_name}:

{pack_citations}

To cite MacSyFinder:

{macsy_cite}
""")
    else:
        print(f"Models '{pack_name}' not found locally.")


def do_check(args: argparse.Namespace) -> None:
    """

    :param args:
    :return:
    """
    # fait a l'instanciation d'un Package

    # si aucune erreur
    # print de la demarche a suivre
    # git tag vers (prendre la version dans les metadata_path)
    # git push --tags remote

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
        description=dedent('''

         *            *               *                   * *       * *      *        *     *
    *           *               *   *   *  *    **                *   * *       *
      **     *    *   *  *     *                    *               *       *
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

    parser.add_argument('-q', '--quiet',
                        action='store_true',
                        help='Give less output.')
    parser.add_argument("-v", "--verbosity",
                        action="count",
                        default=0,
                        help="Give more output.")
    parser.add_argument('-V', '--version',
                        action='version',
                        version='%(prog)s 0.1')

    # -- subparser options -- #

    subparsers = parser.add_subparsers(help=None)

    available_subparser = subparsers.add_parser('available',
                                                help='List Models available on macsy-models')
    available_subparser.add_argument('--org',
                                     default="macsy-models",
                                     help="The name of Model orgagnization"
                                          "(default 'macsy-models'))"
                                     )
    available_subparser.set_defaults(func=do_available)

    download_subparser = subparsers.add_parser('download',
                                               help='Download packages.')

    download_subparser.set_defaults(func=do_download)
    download_subparser.add_argument('-d', '--dest',
                                    help='Download packages into <dir>.')
    download_subparser.add_argument('--org',
                                    default="macsy-models",
                                    help="The name of Model orgagnization"
                                         "(default 'macsy-models'))"
                                    )
    download_subparser.add_argument('package', help='Package name.')
    ###########
    # Install #
    ###########
    install_subparser = subparsers.add_parser('install', help='Install packages.')
    install_subparser.set_defaults(func=do_install)
    install_subparser.add_argument('-f', '--force-reinstall',
                                   help='Reinstall package even if it is already up-to-date.')
    install_subparser.add_argument('--org',
                                   default="macsy-models",
                                   help="The name of Model orgagnization"
                                        "(default 'macsy-models'))"
                                   )
    install_subparser.add_argument('-I', '--ignore-installed',
                                   action='store_true',
                                   default=False,
                                   dest='force',
                                   help='Ignore the installed packages (reinstalling instead).')
    install_subparser.add_argument('-u', '--user',
                                   help='Install to the MacSYFinder user install directory for your platform. '
                                        'Typically ~/.macsyfinder/data')
    install_subparser.add_argument('-U', '--upgrade',
                                   action='store_true',
                                   default=False,
                                   help='Upgrade specified package to the newest available version.')
    install_subparser.add_argument('package',
                                   help='Package name.')
    #############
    # Uninstall #
    #############
    uninstall_subparser = subparsers.add_parser('uninstall',
                                                help='Uninstall packages.')
    uninstall_subparser.set_defaults(func=do_uninstall)
    uninstall_subparser.add_argument('package',
                                     help='Package name.')
    ##########
    # search #
    ##########
    search_subparser = subparsers.add_parser('search',
                                             help='Discover new packages.')
    search_subparser.set_defaults(func=do_search)
    search_subparser.add_argument('--org',
                                  default="macsy-models",
                                  help="The name of Model orgagnization"
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
    ########
    # list #
    ########
    list_subparser = subparsers.add_parser('list',
                                           help='List installed packages.')
    list_subparser.set_defaults(func=do_list)
    list_subparser.add_argument('-o', '--outdated',
                                help='List outdated packages.')
    list_subparser.add_argument('-u', '--uptodate',
                                help='List uptodate packages')
    ########
    # cite #
    ########
    cite_subparser = subparsers.add_parser('cite',
                                           help='How to cite a package.')
    cite_subparser.set_defaults(func=do_cite)
    cite_subparser.add_argument('package',
                                help='Package name.')
    #########
    # check #
    #########
    check_subparser = subparsers.add_parser('check',
                                            help='check if the directory is ready to be publish as data package')
    check_subparser.set_defaults(func=do_check)
    check_subparser.add_argument('dir',
                                 nargs='?',
                                 default=os.getcwd(),
                                 help='the directory to check')
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
    return "macsydata {}".format(func_name)


def main(args=None, loglevel=None) -> None:
    """
    Main entry point.

    :param args: the arguments passed on the command line (before parsing)
    :type args: list
    :rtype: int
    """

    args = sys.argv[1:] if args is None else args
    parser = build_arg_parser()
    parsed_args = parser.parse_args(args)

    defaults = MacsyDefaults()
    config = Config(defaults, argparse.Namespace())

    macsypy.init_logger()
    if not loglevel:
        # logs are specify from args options
        macsypy.logger_set_level(level=config.log_level())
    else:
        # used by unit tests to mute or unmute logs
        macsypy.logger_set_level(level=loglevel)

    logger = logging.getLogger('macsypy.macsydata')

    if 'func' in parsed_args:
        parsed_args.func(parsed_args)
        logger.debug("'{}' command completed successfully.".format(cmd_name(parsed_args)))
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
