#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur (Paris) and CNRS.                         #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


import time
import sys
import os

from ConfigParser import SafeConfigParser, NoSectionError, NoOptionError
try:
    import warnings
except ImportError:
    warnings = None
from distutils import log, dir_util
from distutils.core import setup
from distutils.core import Command
from distutils.dist import Distribution
from distutils.command.build import build
from distutils.command.install import install
from distutils.command.install_data import install_data as _install_data
from distutils.command.build_scripts import build_scripts as _build_scripts
from distutils.command.install_scripts import install_scripts as _install_scripts
from distutils.command.sdist import sdist
from distutils.errors import DistutilsFileError, DistutilsOptionError, DistutilsPlatformError
from distutils.versionpredicate import VersionPredicate
from distutils.util import subst_vars as distutils_subst_vars
from distutils.util import get_platform, change_root, convert_path


class check_and_build( build ):
    
    def run(self):
        chk = True
        for req in require_python:
            chk &= self.check_python(req)
        for req in require_packages:
            chk &= self.check_package(req)
        if not chk: 
            sys.exit(1)
        build.run(self)
        log.info("""
Unit tests are available. It is _highly_ recommended to run tests now.
to run test, run 'python setup.py test -vv'""")

    def check_python(self, req):
        chk = VersionPredicate(req)
        ver = '.'.join([str(v) for v in sys.version_info[:2]])
        if not chk.satisfied_by(ver):
            log.error("Invalid python version, expected {0}".format(req))
            return False
        return True

    def check_package(self, req):
        chk = VersionPredicate(req)
        try:
            mod = __import__(chk.name)
        except:
            log.error("Missing mandatory {0} python module".format(chk.name))
            return False
        for v in ['__version__', 'version']:
            ver = getattr(mod, v, None)
            break
        try:
            if ver and not chk.satisfied_by(ver):
                log.error("Invalid module version, expected {0}".format(req))
                return False
        except:
            pass
        return True

class build_scripts(_build_scripts):

    def finalize_options(self):
        _build_scripts.finalize_options(self)
        self.scripts = self.distribution.scripts 
        self.viewer = self.distribution.viewer

    def run(self):
        self.mkpath(self.build_dir)
        must_clean = []
        for _file in self.distribution.viewer:
            # if the link already exists (build already ran don't try to create it
            if not os.path.exists(_file):
                macsyview_src = os.path.join('macsyview', _file)
                if os.path.exists(macsyview_src):
                    os.symlink(os.path.abspath(os.path.join('macsyview', _file)), _file)
                    must_clean.append(_file)
                    self.scripts += self.viewer
                else:
                    log.warn("macsyview not found")
        _build_scripts.run(self)
        for _link in must_clean:
            os.unlink(_link)


class test(Command):

    description = "run the unit tests against the build library"

    user_options = [('verbosity', 'v', 'verbosity of outputs (cumulative option ex -vv)', 1),
                    ('build-base=', 'b', "base build directory (default: 'build.build-base')"),
                    ('build-lib=', None, "build directory for all modules (default: 'build.build-lib')"),
                    ('plat-name=', 'p', "platform name to build for, if supported (default: {0})".format(get_platform())),
                    ]

    help_options = []

    def initialize_options(self):
        self.verbosity = None
        self.build_base = 'build'
        self.build_lib = None
        self.build_purelib = None
        self.build_platlib = None
        self.plat_name = None
        self.skip_build = 0
        self.warn_dir = 1

    def finalize_options(self):
        # if self.build_lib is None:
        #    self.build_lib = os.path.join(self.build_base, 'lib' )
        if self.verbosity is None:
            self.verbosity = 0
        else:
            self.verbosity = int(self.verbosity)

        if self.plat_name is None:
            self.plat_name = get_platform()
        else:
            # plat-name only supported for windows (other platforms are
            # supported via ./configure flags, if at all).  Avoid misleading
            # other platforms.
            if os.name != 'nt':
                raise DistutilsOptionError(
                            "--plat-name only supported on Windows (try "
                            "using './configure --help' on your platform)")

        plat_specifier = ".{0}-{1}".format(self.plat_name, sys.version[0:3])

        # Make it so Python 2.x and Python 2.x with --with-pydebug don't
        # share the same build directories. Doing so confuses the build
        # process for C modules
        if hasattr(sys, 'gettotalrefcount'):
            plat_specifier += '-pydebug'

        # 'build_purelib' and 'build_platlib' just default to 'lib' and
        # 'lib.<plat>' under the base build directory.  We only use one of
        # them for a given distribution, though --
        if self.build_purelib is None:
            self.build_purelib = os.path.join(self.build_base, 'lib')
        if self.build_platlib is None:
            self.build_platlib = os.path.join(self.build_base,
                                              'lib' + plat_specifier)

        # 'build_lib' is the actual directory that we will use for this
        # particular module distribution -- if user didn't supply it, pick
        # one of 'build_purelib' or 'build_platlib'.
        if self.build_lib is None:
            if os.path.exists(self.build_purelib):
                self.build_lib = self.build_purelib
            elif os.path.exists(self.build_platlib):
                self.build_lib = self.build_platlib


    def run(self):
        """
        """
        if not self.skip_build:
            self.run_command('build')
            # If we built for any other platform, we can't install.
            build_plat = self.distribution.get_command_obj('build').plat_name
            # check warn_dir - it is a clue that the 'install' is happening
            # internally, and not to sys.path, so we don't check the platform
            # matches what we are running.
            if self.warn_dir and build_plat != get_platform():
                raise DistutilsPlatformError("Can't test when "
                                             "cross-compiling")
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from tests.unit import main as unit_main
        if self.build_lib is None:
            if os.path.exists(self.build_purelib):
                self.build_lib = self.build_purelib
            elif os.path.exists(self.build_platlib):
                self.build_lib = self.build_platlib

        result_all_tests = []
        os.environ['MACSY_HOME'] = os.path.dirname(os.path.abspath(__file__))
        kind_of_skipped = {}
        log.info("running unit tests")
        unit_res = unit_main.run(self.build_lib, [], verbosity=self.verbosity)
        for test in unit_res.skipped:
            kind_of_skipped[test[1]] = True
        for skip_reason in kind_of_skipped.keys():
            if skip_reason == "neither makeblast nor formatdb found in PATH":
                msg = """
#####################################################################
neither makeblast nor formatdb found in PATH
You'll have to provide the full path of 
makeblastdb or formatdb in config file or command line to run macsyfinder
#####################################################################"""
            elif skip_reason == "hmmsearch not found in PATH":
                msg = """
########################################################
hmmsearch not found in PATH
You'll have to provide the full path of 
hmmsearch in config file or command line to run macsyfinder
########################################################"""
            else:
                msg = skip_reason
            log.warn(msg)
        result_all_tests.append(unit_res)

        log.info("running functional tests")
        from tests.functional import main as func_main
        func_main = func_main.run(self.build_lib, [], verbosity=self.verbosity)
        result_all_tests.append(func_main)
        sys.path.pop(0)
        if all(result_all_tests):
            sys.exit(0)
        else:
            sys.exit("some tests fails. Run python setup.py test -vv to have more details")



class install_macsyfinder(install):

    # I use record to store all installed files and reuse this record file for uninstall
    # so this option is not available anymore for the users
    for i, opt in enumerate(install.user_options):
        if opt[0] == 'record=':
            install.user_options.pop(i)

    def initialize_options(self):
        install.initialize_options(self)

    def finalize_options(self):
        install.finalize_options(self)
        self.record = self.distribution.uninstall_files
        with open(self.distribution.uninstall_prefix, "w") as _f:
            _f.write("[install]\n")
            _f.write('install_scripts = {}\n'.format(os.path.normpath(self.install_scripts)))
            _f.write('install_lib = {}\n'.format(os.path.normpath(self.install_lib)))


    def run(self):
        def subst_file(_file, vars_2_subst):
            input_file = os.path.join(self.build_lib, _file)
            output_file = input_file + '.tmp'
            subst_vars(input_file, output_file, vars_2_subst)
            os.unlink(input_file)
            self.move_file(output_file, input_file)

        inst = self.distribution.command_options.get('install')
        vars_2_subst = {'PREFIX': inst['prefix'][1] if 'prefix' in inst else '',
                        'PREFIXCONF': os.path.join(get_install_conf_dir(inst), 'macsyfinder'),
                        'PREFIXDATA': os.path.join(get_install_data_dir(inst), 'macsyfinder'),
                        'PREFIXDOC': os.path.join(get_install_doc_dir(inst), 'macsyfinder'),
                        }
        for _file in self.distribution.fix_prefix:
            subst_file(_file, vars_2_subst)
        _file = os.path.join('macsypy', '__init__.py')
        subst_file(_file, {'VERSION': self.distribution.get_version()})
        install.run(self)


class install_data(_install_data):

    user_options = [
        ('install-dir=', 'd',
         "base directory for installing data files "
         "(default: installation base dir)"),
        ('root=', None,
         "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ]

    boolean_options = ['force']


    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        self.install_dir = get_install_data_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                   )
        self.prefix_data = self.install_dir
        self.files_2_install = self.distribution.data_files 
        with open(self.distribution.uninstall_prefix, "a") as _f:
            _f.write('install_data = {}\n'.format(self.install_dir))


    def run(self):
        self.mkpath(self.install_dir)
        for f in self.files_2_install:
            if isinstance(f, str):
                if not os.path.exists(f):
                    log.warn("WARNING the document {} cannot be found, installation skipped".format(f))
                # it's a simple file, so copy it
                f = convert_path(f)
                if self.warn_dir:
                    self.warn("setup script did not provide a directory for "
                              "'{0}' -- installing right in '{1}'".format((f, self.install_dir)))
                (out, _) = self.copy_file(f, self.install_dir)
                self.outfiles.append(out)
            else:
                # it's a tuple with path to install to and a list of path
                dir = convert_path(f[0])
                if not os.path.isabs(dir):
                    dir = os.path.join(self.install_dir, dir)
                elif self.root:
                    dir = change_root(self.root, dir)
                self.mkpath(dir)
                if f[1] == []:
                    # If there are no files listed, the user must be
                    # trying to create an empty directory, so add the
                    # directory to the list of output files.
                    self.outfiles.append(dir)
                else:
                    # Copy files, adding them to the list of output files.
                    for data in f[1]:
                        data = convert_path(data) # return name that will work on the native filesystem
                        if not os.path.exists(data):
                            log.warn("WARNING the document {0} cannot be found, installation skipped".format(data))
                            continue
                        if os.path.isdir(data):
                            out = self.copy_tree(data, dir)
                            self.outfiles.extend(out)
                        else:
                            (out, _) = self.copy_file(data, dir)
                            self.outfiles.append(out)



class install_doc(install_data):

    install.sub_commands += [('install_doc', lambda self: not self.no_doc)]

    description = "installation directory for documentation files"

    setattr(install, 'install_doc', None)
    setattr(install, 'no_doc', None)

    install.user_options.append(('install-doc=', None, description))
    install.user_options.append(('no-doc', None, 'do not install documentation'))

    user_options = [
        ('install-doc=', 'd', "base directory for installing documentation files "
                              "(default: installation base dir share/doc)"),
        ('root=', None, "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ('no-doc', None, 'do not install documentation')
        ]

    boolean_options = ['force']

    def initialize_options(self):
        install_data.initialize_options(self)
        self.install_doc = None
        self.no_doc = None
        self.files_2_install = self.distribution.doc_files 


    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        self.install_dir = get_install_doc_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )
        self.prefix_data = self.install_dir
        self.no_doc = inst.get('no_doc', ('command line', False))[1]
        with open(self.distribution.uninstall_prefix, "a") as _f:
            _f.write('install_doc = {}\n'.format(self.install_dir))

    def run(self):
        install_data.run(self)


class install_conf(install_data):

    install.sub_commands += [('install_conf', lambda self:True)]

    description = "installation directory for configuration files"

    setattr(install, 'install_conf', None)
    install.user_options.append(('install-conf=', None, description)) 

    user_options = [
        ('install-conf=', 'd',
         "base directory for installing configuration files "
         "(default: installation base dir etc)"),
        ('root=', None,
         "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ]

    boolean_options = ['force']

    def initialize_options(self):
        install_data.initialize_options(self)
        self.conf_files = self.distribution.conf_files


    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        self.install_dir = get_install_conf_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                   )
        with open(self.distribution.uninstall_prefix, "a") as _f:
            _f.write('install_conf = {}\n'.format(self.install_dir))

    def run(self):
        self.mkpath(self.install_dir)
        inst = self.distribution.command_options.get('install')
        vars_2_subst = {'PREFIX': inst['prefix'][1] if 'prefix' in inst else '',
                        'PREFIXCONF': os.path.join(get_install_conf_dir(inst), 'macsyfinder'),
                        'PREFIXDATA': os.path.join(get_install_data_dir(inst), 'macsyfinder'),
                        'PREFIXDOC': os.path.join(get_install_doc_dir(inst), 'macsyfinder')
                        }
        for f in self.conf_files:
            if isinstance(f, str):
                # it's a simple file, so copy it
                f = convert_path(f)
                if self.warn_dir:
                    self.warn("setup script did not provide a directory for "
                              "'{0}' -- installing right in '{1}'".format(f, self.install_dir))
                dest = os.path.join(self.install_dir, f + ".new")
                (out, _) = self.copy_file(f, self.install_dir)
                self.outfiles.append(out)
            else:
                # it's a tuple with path to install to and a list of files
                _dir = convert_path(f[0])
                if not os.path.isabs(_dir):
                    _dir = os.path.join(self.install_dir, _dir)
                elif self.root:
                    _dir = change_root(self.root, _dir)
                self.mkpath(_dir)

                if f[1] == []:
                    # If there are no files listed, the user must be
                    # trying to create an empty directory, so add the
                    # directory to the list of output files.
                    self.outfiles.append(_dir)
                else:
                    # Copy files, adding them to the list of output files.
                    for conf in f[1]:
                        conf = convert_path(conf)
                        dest = os.path.join(_dir, os.path.basename(conf) + ".new")
                        (out, _) = self.copy_file(conf, dest)
                        if conf in self.distribution.fix_conf:
                            input_file = out
                            output_file = input_file + '.tmp'
                            subst_vars(input_file, output_file, vars_2_subst)
                            if os.path.exists(input_file):
                                os.unlink(input_file)
                            self.move_file(output_file, input_file)
                            self.outfiles.append(input_file)


class install_viewer(_install_scripts):

    install.sub_commands += [('install_viewer', lambda self: not self.no_viewer)]

    description = 'Install MacSyView, MacSyView allow to view graphically results of MacSyFinder ' \
                  '(need a web browser to work).'

    setattr(install, 'no_viewer', None)
    install.user_options.append(('no-viewer', None, 'do not install MacSyView'))

    user_options = [ ('no-viewer', None, )]

    def initialize_options(self):
        _install_scripts.initialize_options(self)
        self.no_viewer = False

    def finalize_options(self):
        _install_scripts.finalize_options(self)
        inst = self.distribution.command_options.get('install')
        self.no_viewer = inst.get('no_viewer', ('command line', False))[1]


    def run(self):
        def subst(src_file, dest_file, local_vars):
            with open(src_file , 'r') as f_in:
                with open(dest_file, 'w') as f_out:
                    for line in f_in:
                        for var in local_vars:
                            new_line = line.replace('$'+var, local_vars[var])
                        f_out.write(new_line)

        inst = self.distribution.command_options.get('install')
        vars_2_subst = {'MACSYVIEW': os.path.join(get_install_data_dir(inst),
                                                  'macsyfinder',
                                                  'macsyview',
                                                  'app',
                                                  'index.html')}

        for _file in self.distribution.fix_viewer:
            input_file = os.path.join(self.build_dir, _file)
            output_file = input_file + '.tmp'
            subst(input_file, output_file, vars_2_subst)
            os.unlink(input_file)
            self.move_file(output_file, input_file)
        _install_scripts.run(self)



class Uninstall(Command):

    description = "remove installed files"

    user_options = []

    def initialize_options(self):
        self.install_scripts = None
        self.install_lib = None
        self.install_data = None
        self.install_doc = None
        self.install_conf = None


    def finalize_options(self):
        self.parser = SafeConfigParser()
        if not os.path.exists(self.distribution.uninstall_prefix):
            raise DistutilsFileError("Cannot unistall macsyfinder.\n{}: No such file".format(self.distribution.uninstall_prefix))
        used_files = self.parser.read(self.distribution.uninstall_prefix)
        for attr in [ attr for attr in vars(self) if attr.startswith('install_')]:
            try:
                value = self.parser.get('install', attr)
            except(NoSectionError, NoOptionError):
                continue
            setattr(self, attr, value)

    def run(self):
        prefixes = []
        for attr in [ attr for attr in vars(self) if attr.startswith('install_')]:
            prefixes.append(getattr(self, attr))

        def clean_tree(_dir):
            find_prefix = False
            for prefix in prefixes:
                print "== ", prefix, ".find(", _dir, ") = ", prefix.find(_dir)
                if prefix.find(_dir) != -1:
                    find_prefix = True
                    return prefix

            if find_prefix:
                return
            try:
                if not self.dry_run:
                    os.rmdir(_dir)
                log.info("remove dir {}".format(_dir))
            except OSError as err:
                if err.errno == os.errno.ENOTEMPTY:
                    return
                else:
                    self.warn(err)
                    return
            clean_tree(os.path.dirname(_dir))

        try:
            with open(self.distribution.uninstall_files) as record_file:
                for path in record_file:
                    path = os.path.normpath(path.strip())
                    try:
                        if not self.dry_run:
                            os.unlink(path)
                        log.info("remove file {}".format(path))
                    except Exception as err:
                        pass
                    _dir = os.path.dirname(path)
                    clean_tree(_dir)
        except IOError as err:
            msg = "Cannot unistall macsyfinder.\n"
            if err.errno == os.errno.ENOENT:
                msg += "Cannot access \"{0}\": No such file".format(self.distribution.uninstall_files)
            elif err.errno == os.errno.EACCES:
                msg += "Cannot access \"{0}\": Permission denied".format(self.distribution.uninstall_files)
            else:
                msg += str(err)
            raise DistutilsFileError(msg)


###################################
#                                 #
# creation of source distribution #
#                                 #
###################################

class sdist_macsy(sdist):

    user_options = [('macsyview-dir=', None, 'The location of the macsyview project')] + sdist.user_options

    def initialize_options(self):
        sdist.initialize_options(self)
        self.macsyview_dir = None
        self.macsyview_link_name = None

    def run(self):
        inst = self.distribution.command_options.get('sdist')
        self.macsyview_link_name = 'macsyview'
        if not inst.has_key('macsyview_dir'):
            self.must_clean = False
            if not os.path.exists(self.macsyview_link_name):
                sys.exit('Missing mandatory macsyview-dir to specify the location of macsyview project')
        else:
            self.set_up()
        sdist.run(self)
        self.clean_up()

    def update_mktemplate(self):
        with open('MANIFEST.in', 'a') as f:
            f.write('recursive-include {0} *\n'.format(os.path.join(self.macsyview_link_name, 'app')))

    def set_up(self):
        if os.path.exists(self.macsyview_link_name):
            if os.path.realpath(self.macsyview_link_name) == os.path.realpath(self.macsyview_dir):
                self.must_clean = False
                return
            else:
                raise RuntimeError("there already exist a 'macsyview' path "
                                   "in {0} and is different than {1}".format(os.getcwd(),
                                                                             self.macsyview_dir))
        os.symlink(self.macsyview_dir, self.macsyview_link_name)
        self.must_clean = True

    def clean_up(self):
        if self.must_clean:  
            os.unlink(self.macsyview_link_name)


class UsageDistribution(Distribution):

    def __init__(self, attrs=None):
        # It's important to define options before to call __init__
        # otherwise AttributeError: UsageDistribution instance has no attribute 'conf_files'
        self.conf_files = None
        self.doc_files = None
        self.viewer = None
        self.fix_prefix = None
        self.fix_conf = None
        self.fix_viewer = None
        self.uninstall_prefix = os.path.join(os.path.dirname(__file__), "uninstall.cfg")
        self.uninstall_files = os.path.join(os.path.dirname(__file__), "uninstall_files")
        Distribution.__init__(self, attrs = attrs)
        self.common_usage = """\
Common commands: (see '--help-commands' for more)

  setup.py build      will build the package underneath 'build/'
  setup.py test       will run the tests on the newly build library
  setup.py install    will install the package
  setup.py uninstall  will unistall every installed files
"""


def get_install_data_dir(inst):

    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])

    if 'install_data' in inst:
        install_dir = inst['install_data'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share')
    else:
        install_dir = os.path.join('/', 'usr', 'local', 'share')
    return install_dir


def get_install_conf_dir(inst):
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])

    if 'install_conf' in inst:
        install_dir = inst['install_conf'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'etc')
    else:
        install_dir = '/etc'
    return install_dir


def get_install_doc_dir(inst):
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])

    if 'install_doc' in inst:
        install_dir = inst['install_doc'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share', 'doc')
    else:
        install_dir = os.path.join('/', 'usr', 'local', 'share', 'doc')
    return install_dir


def subst_vars(src, dst, vars):
    try:
        src_file = open(src, "r")
    except os.error as err:
        raise DistutilsFileError, "could not open '{0}': {1)".format(src, err)
    try:
        dest_file = open(dst, "w")
    except os.error as err:
        raise DistutilsFileError, "could not create '{0}': {1}".format(dst, err)
    with src_file:
        with dest_file:
            for line in src_file:
                new_line = distutils_subst_vars(line, vars)
                dest_file.write(new_line)


require_python = ['python (>=2.7, <3.0)']
require_packages = []


setup(name='macsyfinder',
      version='1.0.4',
      description="""MacSyFinder: Detection of macromolecular systems 
in protein datasets using systems modelling and similarity search""",
      author="Sophie Abby, Bertrand Néron",
      author_email="sabby@pasteur.fr, bneron@pasteur.fr",
      classifiers=[
                  'Operating System :: POSIX',
                  'Programming Language :: Python :: 2.7 :: Only',
                  'Topic :: Bioinformatics',
                 ],
      packages=['macsypy'],
      scripts=['bin/macsyfinder'],
      viewer=['bin/macsyview'],
      # (dataprefix +'where to put the data in the install, [where to find the data in the tar ball]
      data_files=[('macsyfinder/DEF', ['data/DEF/']),
                  ('macsyfinder/profiles', ['data/profiles/']),
                  ('macsyfinder/sequence_data', ['data/sequence_data/']),
                  ('macsyfinder/macsyview/app', ['macsyview/app/'])
                 ],
      conf_files=[('macsyfinder', ['etc/macsyfinder.conf'])],
      doc_files=[('macsyfinder/html', ['doc/_build/html/']),
                 ('macsyfinder/pdf', ['doc/_build/latex/Macsyfinder.pdf']),
                ],
      # file where some variable must be fix by install_conf
      fix_conf=['etc/macsyfinder.conf'],
      # file where some variable must be fix by macsyfinder_install
      fix_prefix=['macsypy/config.py', 'macsypy/registries.py'],
      fix_viewer=['macsyview'],
      cmdclass={'build': check_and_build,
                'build_scripts': build_scripts,
                'test': test,
                'install': install_macsyfinder,
                'install_viewer': install_viewer,
                'install_data': install_data,
                'install_conf': install_conf,
                'install_doc': install_doc,
                'uninstall': Uninstall,
                'sdist': sdist_macsy,
                },
      distclass=UsageDistribution
      )

