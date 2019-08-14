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
import shutil
from textwrap import dedent
from typing import List, Tuple, Dict
import pathlib

import colorlog
from packaging import requirements
from packaging import specifiers

import macsypy
from macsypy.config import MacsyDefaults, Config
from macsypy.registries import ModelRegistry, scan_models_dir
from macsypy.package import RemoteModelIndex, Package


def init_logger(level='INFO', out=True):
    logger = colorlog.getLogger('macsydata')
    logging = colorlog.logging.logging
    handlers = []
    if out:
        stdout_handler = colorlog.StreamHandler(sys.stderr)
        stdout_formatter = colorlog.ColoredFormatter("%(log_color)s%(message)s",
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

    logger.setLevel(getattr(logging, level))
    return logger


_log = init_logger()

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


def _search_in_pack_name(pattern: str, remote: RemoteModelIndex, packages: List[str],
                         match_case: bool = False) -> List[Tuple[str, str, Dict]]:
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
            _log.info(f"Downloading {pack_name} {vers}")
            arch_path = remote.download(pack_name, vers, dest=args.dest)
            _log.info(f"Successfully downloaded packaging {pack_name} in {arch_path}")
            return arch_path
        else:
            _log.error(f"No version that satisfy requirements '{specifier}' for '{pack_name}'.")
            _log.warning(f"Available versions: {','.join(all_versions)}")


def _find_all_packages() -> ModelRegistry:
    """
    :return: all madels installed
    """
    defaults = MacsyDefaults()
    config = Config(defaults, argparse.Namespace())
    system_model_dir = config.models_dir()
    user_model_dir = os.path.join(os.path.expanduser('~'), '.macsyfinder', 'data')
    model_dirs = (system_model_dir, user_model_dir) if os.path.exists(user_model_dir) else (system_model_dir,)
    registry = ModelRegistry()
    for model_dir in model_dirs:
        for model_loc in scan_models_dir(model_dir, profile_suffix=config.profile_suffix):
            registry.add(model_loc)
    return registry


def _find_local_package(pack_name):
    """

    :param pack_name: the name of the family model to search
    :return: The model location corresponding to the `pack_name`
    :rtype: :class:`macsypy.registries.ModelLocation` object
    """
    registry = _find_all_packages()
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
    local_spec = False
    if local_pack:
        pack = Package(local_pack.path)
        local_vers = pack.metadata['vers']

    remote = RemoteModelIndex(org=args.org)
    all_versions = remote.list_package_vers(pack_name)
    if not specifier and local_pack:
        specifier = specifiers.SpecifierSet(f">{local_vers}")
        local_spec = True
    compatible_version = list(specifier.filter(all_versions))
    if local_pack and compatible_version and not args.upgrade:
        _log.warning(f"{pack_name} is already installed but {compatible_version[0]} is available.\n"
                     f"To install it please run 'macsydata install --upgrade {pack_name}'")
        return None
    elif local_pack and not compatible_version:
        if args.force:
            if local_vers in all_versions:
                target_vers = local_vers
            elif all_versions:
                target_vers = all_versions[0]
            else:
                _log.warning(f"Cannot find version for models '{pack_name}'.")
                return None
        elif local_spec:
            _log.warning(f"Requirement already satisfied: {pack_name}{specifier} in {pack.path}.\n"
                         f"To force installation use option -f --force-reinstall.")
            return None
        else:
            _log.error(f"Could not find version that satisfied '{pack_name}{specifier}'")
            return None
    elif not local_pack and not compatible_version:
        _log.error(f"Could not find version that satisfied '{pack_name}{specifier}'")
        return None
    else:
        target_vers = compatible_version[0]

    _log.info(f"Downloading {pack_name} ({target_vers}).")
    tmp_arch = remote.download(pack_name, target_vers)
    _log.info(f"Extracting {pack_name} ({target_vers}).")
    cached_pack = remote.unarchive_package(tmp_arch)
    if args.user:
        dest = os.path.realpath(os.path.join(os.path.expanduser('~'), '.macsyfinder', 'data'))
        if os.path.exists(dest) and not os.path.isdir(dest):
            raise RuntimeError("'{}' already exist and is not a directory.")
        elif not os.path.exists(dest):
            os.makedirs(dest)
    else:
        defaults = MacsyDefaults()
        config = Config(defaults, argparse.Namespace())
        dest = config.models_dir()
    if local_pack:
        old_pack_path = f"{local_pack.path}.old"
        shutil.move(local_pack.path, old_pack_path)
    _log.info(f"Installing {pack_name} ({target_vers}) in {dest}")
    shutil.move(cached_pack, dest)
    if local_pack:
        _log.info("Cleaning.")
        shutil.rmtree(old_pack_path)
        shutil.rmtree(pathlib.Path(cached_pack).parent)
    _log.info(f"The models {pack_name} ({target_vers}) have been installed successfully.")


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
        pack = Package(local_pack.path)
        shutil.rmtree(pack.path)
        _log.info(f"models '{pack_name}' in {pack.path} uninstalled.")
    else:
        _log.error(f"Models '{pack_name}' not found locally.")


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
        pack = Package(local_pack.path)
        print(pack.info())
    else:
        _log.error(f"Models '{pack_name}' not found locally.")


def do_list(args: argparse.Namespace) -> None:
    """
    List installed models.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    registry = _find_all_packages()
    for model_loc in registry.models():
        try:
            pack = Package(model_loc.path)
            pack_vers = pack.metadata['vers']
            if args.outdated or args.uptodate:
                remote = RemoteModelIndex(org=args.org)
                all_versions = remote.list_package_vers(pack.name)
                specifier = specifiers.SpecifierSet(f">{pack_vers}")
                update_vers = list(specifier.filter(all_versions))
                if args.outdated and update_vers:
                    print(f"{model_loc.name}-{update_vers}) [{pack_vers}]")
                if args.uptodate and not update_vers:
                    print(f"{model_loc.name}-{pack_vers}")
            else:
                print(f"{model_loc.name}-{pack_vers}")
        except Exception as err:
            if args.verbose > 1:
                _log.warning(str(err))


def do_freeze(args: argparse.Namespace) -> None:
    """
    """
    registry = _find_all_packages()
    for model_loc in sorted(registry.models(), key=lambda ml: ml.name.lower()):
        try:
            pack = Package(model_loc.path)
            pack_vers = pack.metadata['vers']
            print(f"{model_loc.name}=={pack_vers}")
        except:
            pass


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
        pack = Package(local_pack.path)
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
        _log.error(f"Models '{pack_name}' not found locally.")


def do_check(args: argparse.Namespace) -> None:
    """

    :param args:
    :return:
    """
    pack = Package(args.path)
    errors, warnings = pack.check()
    if errors:
        for error in errors:
            _log.error(error)
        _log.error("Please fix issues above, before publishing these models.")
        sys.exit(5)
    if warnings:
        for warning in warnings:
            _log.warning(warning)
        _log.warning("It is better, if you fix warnings above, before to publish these models.")

    if not warnings:
        _log.info("If everyone were like you, I'd be out of business")
        _log.info("To push the models in organization:")
        if os.path.realpath(os.getcwd()) != pack.path:
            # I use level 25 just to remove color
            _log.log(25, f"\tcd {pack.path}")
        if not os.path.exists(os.path.join(pack.path, '.git')):
            _log.info("Transform the models into a git repository ")
            _log.log(25, "\tgit init .")
            _log.log(25, "\tgit add .")
            _log.log(25, "\tgit commit -m 'initial commit'")
            _log.info("add a remote repository to host the models")
            _log.info("for instance if you want to add the models to 'macsy-models'")
            _log.log(25, "\tgit remote add origin https://github.com/macsy-models/")

        _log.log(25, f"\tgit tag {pack.metadata['vers']}")
        _log.log(25, f"\tgit push --tags")

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
    parser.add_argument("-v", "--verbose",
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
                                    help='Download packages into <dir>.')
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
    install_subparser.add_argument('-I', '--ignore-installed',
                                   action='store_true',
                                   default=False,
                                   dest='force',
                                   help='Ignore the installed packages (reinstalling instead).')
    install_subparser.add_argument('-u', '--user',
                                   action='store_true',
                                   default=False,
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
    ##########
    # freeze #
    ##########
    freeze_subparser = subparsers.add_parser('freeze',
                                             help='List installed models in requirements format.')
    freeze_subparser.set_defaults(func=do_freeze)
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
    check_subparser.add_argument('path',
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


def main(args=None) -> None:
    """
    Main entry point.

    :param args: the arguments passed on the command line (before parsing)
    :type args: list
    :rtype: int
    """
    args = sys.argv[1:] if args is None else args
    parser = build_arg_parser()
    parsed_args = parser.parse_args(args)

    macsypy.init_logger()

    if 'func' in parsed_args:
        parsed_args.func(parsed_args)
        _log.debug("'{}' command completed successfully.".format(cmd_name(parsed_args)))
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
