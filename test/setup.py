#!/usr/bin/env python
# -*- coding: utf-8 -*-
from distutils import log
from distutils.core import setup
#from distutils.command.install_egg_info import install_egg_info
from distutils.command.build import build
from distutils.command.install import install
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
        build.run( self )

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

class install_txsscan(install):

    def run(self):
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
require_packages = ['bsddb3']
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
                  'install' : install_txsscan
                 }
      )
