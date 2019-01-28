#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur, Paris.                                   #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


import sys
import os
import argparse
import logging
from textwrap import dedent


def do_download(args):
    """
    Download tarball from remote models repository.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    raise Exception('Not implemented')


def do_install(args):
    """
    Install new models in macsyfinder local models repository.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    raise Exception('Not implemented')


def do_uninstall(args):
    """
    Remove models from macsyfinder local models repository.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    raise Exception('Not implemented')


def do_search(args):
    """
    Models discovery.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    raise Exception('Not implemented')


def do_show(args):
    """
    Show model details.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    raise Exception('Not implemented')


def do_list(args):
    """
    List installed models.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    raise Exception('Not implemented')


def do_cite(args):
    """
    How to cite macsyfinder.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    raise Exception('Not implemented')


def do_help(args):
    """
    Positional help command (alias for '--help' option).

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    if args.topic is None:
        parser.print_help()
    else:
        if args.topic in subparsers.choices:
            subparsers.choices[args.topic].print_help()
        else:
            print('Help topic not found (%s)'%args.topic, file=sys.stderr)


def build_arg_parser():
    """
    Build argument parser.

    :rtype: :class:`argparse.ArgumentParser` object
    """
    global parser, subparsers # to be accessible from 'do_help' func

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

    parser.add_argument('-q','--quiet', action='store_true', help='Give less output.')
    parser.add_argument("-v", "--verbosity", action="count", dest="verbosity", default=0, help="Give more output.")
    parser.add_argument('-V', '--version', action='version', version='%(prog)s 0.1')

    # -- subparser options -- #

    subparsers = parser.add_subparsers(help=None)

    subparser = subparsers.add_parser('download', help='Download packages.')
    subparser.set_defaults(func=do_download)
    subparser.add_argument('-d', '--dest', help='Download packages into <dir>.')
    subparser.add_argument('-i', '--index-url', help='Base URL of macsyfinder package index (default https://github.com/gem/macsyfinder).')
    subparser.add_argument('package', help='Package name.')

    subparser = subparsers.add_parser('install', help='Install packages.')
    subparser.set_defaults(func=do_install)
    subparser.add_argument('-f', '--force-reinstall', help='Reinstall package even if it is already up-to-date.')
    subparser.add_argument('-i', '--index-url', help='')
    subparser.add_argument('-I', '--ignore-installed', help='Ignore the installed packages (reinstalling instead).')
    subparser.add_argument('-t', '--target', help='Install packages into <dir>.')
    subparser.add_argument('-u', '--user', help='Install to the Python user install directory for your platform. Typically ~/.local/.')
    subparser.add_argument('-U', '--upgrade', help='Upgrade specified package to the newest available version.')
    subparser.add_argument('package', help='Package name.')

    subparser = subparsers.add_parser('uninstall', help='Uninstall packages.')
    subparser.set_defaults(func=do_uninstall)
    subparser.add_argument('-y', '--yes', help="Don't ask for confirmation of uninstall deletions.")
    subparser.add_argument('package', help='Package name.')

    subparser = subparsers.add_parser('search', help='Discover new packages.')
    subparser.set_defaults(func=do_search)
    subparser.add_argument('-i', '--index-url', help='')
    subparser.add_argument('pattern', help='Searches for packages matching the pattern.')

    subparser = subparsers.add_parser('show', help='Show information about packages.')
    subparser.set_defaults(func=do_show)
    subparser.add_argument('-f', '--files', help='Show the full list of installed files for each package.')
    subparser.add_argument('-i', '--index-url', help='')
    subparser.add_argument('package', help='Package name.')

    subparser = subparsers.add_parser('list', help='List installed packages.')
    subparser.set_defaults(func=do_list)
    subparser.add_argument('-o', '--outdated', help='List outdated packages.')
    subparser.add_argument('-u', '--uptodate', help='List uptodate packages')

    subparser = subparsers.add_parser('cite', help='How to cite a package.')
    subparser.set_defaults(func=do_cite)
    subparser.add_argument('package', help='Package name.')

    subparser = subparsers.add_parser('help', help='Show help.')
    subparser.set_defaults(func=do_help)
    subparser.add_argument('topic', nargs='?')

    return parser


def log_init(args):
    """
    Log initialization.

    :param args: the arguments passed on the command line
    :type args: :class:`argparse.Namespace` object
    :rtype: None
    """
    global logger # to be accessible from anywhere in the script

    FORMAT_SH = '%(message)s (L %(lineno)d)'
    FORMAT_FH = '%(asctime)-15s : %(levelname)-8s : L %(lineno)d : %(message)s'

    if args.verbosity == 0:
        log_level = logging.ERROR
    elif args.verbosity == 1:
        log_level = logging.WARNING
    elif args.verbosity == 2:
        log_level = logging.INFO
    elif args.verbosity == 3:
        log_level = logging.DEBUG
    else:
        log_level = logging.NOTSET

    logger = logging.getLogger()
    logger.setLevel(log_level)

    sh = logging.StreamHandler(sys.stderr)
    sh.setLevel(log_level)
    sh.setFormatter(logging.Formatter(FORMAT_SH))

    fh = logging.FileHandler('macsydata.log')
    fh.setLevel(log_level)
    fh.setFormatter(logging.Formatter(FORMAT_FH))

    logger.addHandler(fh)
    logger.addHandler(sh)


def cmd_name(args):
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
    func_name = args.func.__name__.replace('do_','')
    return "macsydata {}".format(func_name)


def main(args=None):
    """
    Main entry point.

    :param args: the arguments passed on the command line (before parsing)
    :type args: list
    :rtype: int
    """

    try:
        args = sys.argv[1:] if args is None else args
        parser = build_arg_parser()
        parsed_args = parser.parse_args(args)

        log_init(parsed_args)

        if 'func' in parsed_args:
            parsed_args.func(parsed_args)
            logger.debug("'{}' command completed successsfully.".format(cmd_name(parsed_args)))
        else:
            parser.print_help()

        return 0
    except Exception as e:
        if 'logger' in globals():
            logger.error("'{}' command failed.".format(cmd_name(parsed_args)))
            logger.error("{}".format(str(e)))
        return 1


if __name__ == "__main__":
    sys.exit(main())
