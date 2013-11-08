#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import sys
import os

from distutils import log
from distutils.core import setup
from distutils.dist import Distribution
from distutils.core import Command
#from distutils.command.install_egg_info import install_egg_info
from distutils.command.build import build
from distutils.command.install import install
from distutils.command.install_data import install_data as _install_data
from distutils.util import get_platform
from distutils.errors import DistutilsOptionError, DistutilsPlatformError
from distutils.versionpredicate import VersionPredicate
from distutils.errors import DistutilsFileError
from distutils.util import subst_vars as distutils_subst_vars
from distutils.util import change_root, convert_path


class check_and_build( build ):
    def run(self):
        chk = True
        for req in require_python:
            chk &= self.check_python(req)
        for req in require_packages:
            chk &= self.check_package(req)
        if not chk: 
            sys.exit(1)
        test_res_path = os.path.join(self.build_base, ".tests_results")
        try:
            os.unlink(test_res_path)
        except OSError, err:
            pass
        build.run(self)
        print """
Unit tests are available. It is _highly_ recommended to run tests now.
to run test, run 'python setup.py test -vv'"""

    def check_python(self, req):
        chk = VersionPredicate(req)
        ver = '.'.join([str(v) for v in sys.version_info[:2]])
        if not chk.satisfied_by(ver):
            print >> sys.stderr, "Invalid python version, expected %s" % req
            return False
        return True

    def check_package(self, req):
        chk = VersionPredicate(req)
        try:
            mod = __import__(chk.name)
        except:
            print >> sys.stderr, "Missing mandatory %s python module" % chk.name
            return False
        for v in [ '__version__', 'version' ]:
            ver = getattr(mod, v, None)
            break
        try:
            if ver and not chk.satisfied_by(ver):
                print >> sys.stderr, "Invalid module version, expected %s" % req
                return False
        except:
            pass
        return True


class test(Command):

    description = "run the unit tests against the build library"

    user_options = [('verbosity' , 'v' , 'verbosity of outputs (cumulative option ex -vv)', 1),
                    ('build-base=', 'b', "base build directory (default: 'build.build-base')"),
                    ('build-lib=', None, "build directory for all modules (default: 'build.build-lib')"),
                    ('plat-name=', 'p', "platform name to build for, if supported (default: %s)" % get_platform()),
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
        #if self.build_lib is None:
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

        plat_specifier = ".%s-%s" % (self.plat_name, sys.version[0:3])

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
        from test import main
        if self.build_lib is None:
            if os.path.exists(self.build_purelib):
                self.build_lib = self.build_purelib
            elif os.path.exists(self.build_platlib):
                self.build_lib = self.build_platlib

        print "running test"
        os.environ['TXSSCAN_HOME'] = os.path.dirname(os.path.abspath(__file__))
        test_res = main.run(self.build_lib, [], verbosity = self.verbosity)
        res_path = os.path.join(self.build_base, ".tests_results")
        with open(res_path, 'w') as _file:
            print >> _file, int(test_res.wasSuccessful())
        if not test_res.wasSuccessful():
            sys.exit("some tests fails. Run python setup.py test -vv to have more details")



class install_txsscan(install):

    def run(self):
        test_res_path = os.path.join(self.build_base, ".tests_results")
        test_res = 0 # test fails
        if os.path.exists(test_res_path):
            with open(test_res_path) as _file:
                test_res = int(_file.readline().strip())
                if not test_res:
                    msg = "Unit tests failed. It is _highly_ recommended to fix the problem, before installing"
        else:
            msg = """Unit tests are available. It is _highly_ recommended to run tests now, before installing
to run test, run 'python setup.py test'"""

        if not test_res: #test_res = 0 => test fails or test not ran
            test_OK = raw_input( "{}\nAre you sure you want to install anyway (y/N) ?".format(msg))
            if test_OK.lower() in ('y', 'yes'):
                test_OK = True
            else:
                test_OK = False
        else:
            test_OK = True
        if test_OK:
            inst = self.distribution.command_options.get('install')
            vars_2_subst = {'PREFIX': inst.get('prefix', ''),
                            'PREFIXCONF' : os.path.join(get_install_conf_dir(inst), 'txsscan'),
                            'PREFIXDATA' : os.path.join(get_install_data_dir(inst), 'txsscan'),
                            'PREFIXDOC'  : os.path.join(get_install_doc_dir(inst), 'txsscan')
                            }
            for _file in fix_prefix:
                input_file = os.path.join(self.build_lib, _file)
                output_file =  input_file + '.tmp'
                subst_vars(input_file, output_file, vars_2_subst)
                os.unlink(input_file)
                self.move_file(output_file, input_file)
            install.run(self)


class install_data(_install_data):

    #install.sub_commands += [('install_data', lambda self:True)]

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

    def run(self):
        self.mkpath(self.install_dir)
        for f in self.files_2_install:
            if isinstance(f, str):
                # it's a simple file, so copy it
                f = convert_path(f)
                if self.warn_dir:
                    self.warn("setup script did not provide a directory for "
                              "'%s' -- installing right in '%s'" %
                              (f, self.install_dir))
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
                        data = convert_path(data)#return name that will work on the native filesystem
                        if os.path.isdir(data):
                            out = self.copy_tree(data, dir)
                            self.outfiles.extend(out)
                        else:
                            (out, _) = self.copy_file(data, dir)
                            self.outfiles.append(out)



class install_doc(install_data):

    install.sub_commands += [('install_doc', lambda self:True)]

    description = "installation directory for documentation files"

    setattr(install, 'install_doc', None)

    install.user_options.append(('install-doc=', None, description))
    install.user_options.append(('no-doc', None, 'do not install documentation'))

    user_options = [
        ('install-doc=', 'd', "base directory for installing documentation files " "(default: installation base dir share/doc)"),
        ('root=', None, "install everything relative to this alternate root directory"),
        ('force', 'f', "force installation (overwrite existing files)"),
        ('no-doc', None, 'do not install documentation')
        ]

    boolean_options = ['force']

    def initialize_options(self):
        install_data.initialize_options(self)
        self.install_doc = None
        self.no_doc = False
        self.files_2_install = doc_files #as defined at the global level of this file
        
    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        self.install_dir = get_install_doc_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )
        self.prefix_data = self.install_dir
        



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
        self.conf_files = conf_files #as defined at the top of this file


    def finalize_options(self):
        inst = self.distribution.command_options.get('install')
        self.install_dir = get_install_conf_dir(inst)
        self.set_undefined_options('install',
                                   ('root', 'root'),
                                   ('force', 'force'),
                                  )

    def run(self):
        self.mkpath(self.install_dir)
        inst = self.distribution.command_options.get('install')
        vars_2_subst = {'PREFIX': inst['prefix'][1] if 'prefix' in inst else '',
                        'PREFIXCONF' : os.path.join(get_install_conf_dir(inst), 'txsscan'),
                        'PREFIXDATA' : os.path.join(get_install_data_dir(inst), 'txsscan'),
                        'PREFIXDOC'  : os.path.join(get_install_doc_dir(inst), 'txsscan')
                        }
        for f in self.conf_files:
            if isinstance(f, str):
                # it's a simple file, so copy it
                f = convert_path(f)
                if self.warn_dir:
                    self.warn("setup script did not provide a directory for "
                              "'%s' -- installing right in '%s'" %
                              (f, self.install_dir))
                (out, _) = self.copy_file(f, self.install_dir)
                self.outfiles.append(out)
            else:
                # it's a tuple with path to install to and a list of files
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
                    for conf in f[1]:
                        conf = convert_path(conf)
                        (out, _) = self.copy_file(conf, dir)
                        if conf in fix_conf:
                            input_file = out
                            output_file =  input_file + '.tmp'
                            subst_vars(input_file, output_file, vars_2_subst)
                            os.unlink(input_file)
                            self.move_file(output_file, input_file)
                        self.outfiles.append(out)





class UsageDistribution(Distribution):

    def __init__(self,  attrs=None):
        Distribution.__init__(self, attrs = attrs)
        self.common_usage = """\
Common commands: (see '--help-commands' for more)

  setup.py build      will build the package underneath 'build/'
  setup.py test       will run the tests on the newly build library
  setup.py install    will install the package
"""

def get_install_data_dir(inst):
    
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])
    
    if 'install_data' in inst:
        install_dir = inst['install_data'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share')
    else:
        install_dir = os.path.join('/', 'usr', 'share')
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
        install_dir = os.path.join(inst['prefix'][1], 'share', 'doc' )
    else:
        install_dir = os.path.join('/', 'usr', 'share', 'doc')
    return install_dir

        
def subst_vars(src, dst, vars):
    try:
        src_file = open(src, "r")
    except os.error, err:
        raise DistutilsFileError, "could not open '%s': %s" % (src, err)
    try:
        dest_file = open(dst, "w")
    except os.error, err:
        raise DistutilsFileError, "could not create '%s': %s" % (dst, err)
    with src_file:
        with dest_file:
            for line in src_file:
                new_line = distutils_subst_vars(line, vars)
                dest_file.write(new_line)

require_python = [ 'python (>=2.7, <3.0)' ]
require_packages = []

#I cannot succeed to inject conf_file in a distribution
#so i put it at the top level :-(
conf_files = [('txsscan', ['etc/txsscan.conf'])]
doc_files = [('txsscan', ['doc/_build/html/'])]

#file where some variable must be fix by install_conf
fix_conf = ['etc/txsscan.conf']

#file where some variable must be fix by txsscan_install
fix_prefix = ['txsscanlib/config.py', 'txsscanlib/registries.py']


setup(name        = 'txsscan',
      version     =  time.strftime("%Y%m%d"),
      description  = """ """,
      classifiers = [
                     'Operating System :: POSIX' ,
                     'Programming Language :: Python' ,
                     'Topic :: Bioinformatics' ,
                    ] ,
      packages    = ['txsscanlib'],
      scripts     = [ 'bin/txsscan' ] ,
      data_files = [('txsscan/DEF', ['data/DEF/']),
                    ('txsscan/profiles', ['data/profiles/'])
                    ],
      cmdclass= { 'build' : check_and_build ,
                  'test': test,
                  'install' : install_txsscan,
                  'install_data' : install_data,
                  'install_conf' : install_conf,
                  'install_doc'  : install_doc
                 }
      )

