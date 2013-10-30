#!/usr/bin/env python
# -*- coding: utf-8 -*-
from distutils import log
from distutils.core import setup
from distutils.core import Command
#from distutils.command.install_egg_info import install_egg_info
from distutils.command.build import build
from distutils.command.install import install
from distutils.util import get_platform
from distutils.errors import DistutilsOptionError, DistutilsPlatformError
from distutils.versionpredicate import VersionPredicate
from distutils.errors import DistutilsFileError
from distutils.util import subst_vars as distutils_subst_vars

import time
import sys
import os

class check_and_build( build ):
    def run(self):
        chk = True
        for req in require_python:
            chk &= self.check_python(req)
        for req in require_packages:
            chk &= self.check_package(req)
        if not chk: 
            sys.exit(1)
        test_res_path = os.path.join(self.build_lib, ".tests_results")
        try:
            os.unlink(test_res_path)
        except OSError:
            pass
        build.run(self)
        print """
Unit tests are available. It is _highly_ recommended to run tests now.
to run test, run 'python setup.py test'"""

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
        sys.path.insert(0, os.path.join(os.getcwd(), 'test'))
        import main
        if self.build_lib is None:
            if os.path.exists(self.build_purelib):
                self.build_lib = self.build_purelib
            elif os.path.exists(self.build_platlib):
                self.build_lib = self.build_platlib

        print "running test"
        test_res = main.run(self.build_lib, [], verbosity = self.verbosity)
        res_path = os.path.join(self.build_lib, ".tests_results")

        with open(res_path, 'w') as _file:
            print >> _file, int(test_res.wasSuccessful())
        if not test_res.wasSuccessful():
            sys.exit("some tests fails. Run python setup.py test -vv to have more details")



class install_txsscan(install):

    def run(self):
        test_res_path = os.path.join(self.build_lib, ".tests_results")
        test_res = 0
        if os.path.exists(test_res_path):
            with open(test_res_path) as _file:
                test_res = int(_file.readline().strip())
                if not test_res:
                    msg = "Unit tests failed. It is _highly_ recommended to fix the problem, before installing"
        else:
            msg = """Unit tests are available. It is _highly_ recommended to run tests now, before installing
to run test, run 'python setup.py test'"""

        if not test_res: #test_res = 0 => test fails ore test not ran
            test_OK = raw_input( "{}\nAre you sure you want to install anyway (y/N) ?".format(msg))
            if test_OK.lower() in ('y', 'yes'):
                test_OK = True
            else:
                test_OK = False
        else:
            #test_res = 1
            test_OK = True
        if test_OK:
            for _file in fix_prefix:
                input_file = os.path.join(self.build_lib, _file)
                output_file =  input_file + '.tmp'
                subst_vars( input_file, output_file, {'PREFIX': self.prefix})
                os.unlink(input_file)
                self.move_file(output_file, input_file) 
                install.run(self)


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
fix_prefix = ["txsscanlib/config.py"]

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
      data_files=[('etc/txsscan', ['etc/txsscan.conf'])],
      cmdclass= { 'build' : check_and_build ,
                  'test': test,
                  'install' : install_txsscan
                 }
      )
