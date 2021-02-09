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
import sysconfig
import warnings

from distutils.errors import DistutilsFileError, DistutilsSetupError
from distutils.util import subst_vars as distutils_subst_vars

from setuptools import setup, find_packages
from setuptools.dist import Distribution
from setuptools.command.install_lib import install_lib as _install_lib


import sys
if 'setup.py' in sys.argv:
    INSTALLER = 'setuptools'
else:
    INSTALLER = 'pip'

from macsypy import __version__ as mf_version


class install_lib(_install_lib):

    def finalize_options(self):
        inst = self.distribution.command_options.get('install', {})
        _install_lib.finalize_options(self)

    def run(self):
        if INSTALLER == 'pip':
            script_tmp_dir = self.build_dir
        else:  # setuptools
            raise DistutilsSetupError("'setuptools' is not supported. Please use 'pip' instead.")

        def subst_file(_file, vars_2_subst):
            input_file = os.path.join(script_tmp_dir, _file)
            output_file = input_file + '.tmp'
            subst_vars(input_file, output_file, vars_2_subst)
            os.unlink(input_file)
            self.move_file(output_file, input_file)

        inst = self.distribution.command_options.get('install', {})
        _file = os.path.join('macsypy', '__init__.py')
        subst_file(_file, {'MACSYDATA': os.path.join(get_install_data_dir(inst), 'macsyfinder'),
                           'MACSYCONF': os.path.join(get_install_conf_dir(inst), 'macsyfinder')})

        _install_lib.run(self)


class UsageDistribution(Distribution):

    def __init__(self, attrs=None):
        # It's important to define options before to call __init__
        # otherwise AttributeError: UsageDistribution instance has no attribute 'conf_files'
        self.fix_prefix = None
        Distribution.__init__(self, attrs=attrs)
        self.common_usage = """\
Common commands: (see '--help-commands' for more)

  setup.py build      will build the package underneath 'build/'
  setup.py install    will install the package
  setup.py test       run tests after in-place build
"""


def get_install_data_dir(inst):
    """
    :param inst: installation option
    :type inst: dict
    :return: the prefix where to install data
    :rtype: string
    """
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])
    elif 'user' in inst:
        import site
        inst['prefix'] = ('command line', site.USER_BASE)
    elif 'root' in inst:
        inst['prefix'] = ('command line',
                          os.path.join(inst['root'][1],
                                       sysconfig.get_path('data').strip(os.path.sep)
                                       )
                          )

    if 'install_data' in inst:
        install_dir = inst['install_data'][1]
    elif 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'share')
    else:
        try:
            from pip.locations import distutils_scheme
        except ImportError:
            # from pip >=10 distutils_scheme has move into _internal package.
            # It's ugly but I haven't other choice
            # because the asshole Debian/Ubuntu patch pip to install in /usr/local
            # but not python. So when use sysconfig.get_paths['data'] it return '/usr'
            from pip._internal.locations import distutils_scheme
        install_dir = os.path.join(distutils_scheme('')['data'], 'share')
    return install_dir


def get_install_conf_dir(inst):
    """
    :param inst: installation option
    :type inst: dict
    :return: the prefix where to install data
    :rtype: string
    """
    if 'VIRTUAL_ENV' in os.environ:
        inst['prefix'] = ('environment', os.environ['VIRTUAL_ENV'])
    elif 'user' in inst:
        import site
        inst['prefix'] = ('command line', site.USER_BASE)
    elif 'root' in inst:
        inst['prefix'] = ('command line',
                          os.path.join(inst['root'][1], 'etc')
                          )

    if 'prefix' in inst:
        install_dir = os.path.join(inst['prefix'][1], 'etc')
    else:
        install_dir = '/etc'
    return install_dir


def subst_vars(src, dst, vars):
    """
    substitute variables (string starting with $) in file
    :param src: the file containing variable to substitute
    :type src: string
    :param dst: the destination file
    :type dst: string
    :param vars: the variables to substitute in dict key are variable name
    :type vars: dict
    """
    try:
        src_file = open(src, "r")
    except os.error as err:
        raise DistutilsFileError("could not open '{0}': {1}".format(src, err))
    try:
        dest_file = open(dst, "w")
    except os.error as err:
        raise DistutilsFileError("could not create '{0}': {1}".format(dst, err))
    with src_file:
        with dest_file:
            try:
                file = ''
                for line in src_file:
                    file += line
                    new_line = distutils_subst_vars(line, vars)
                    dest_file.write(new_line)
            except UnicodeDecodeError as err:
                import sys
                print("####################### file raising a error#######################3", file=sys.stderr)
                print(file, file=sys.stderr)
                print("###################### end of debug ###############################", file=sys.stderr)
                raise RuntimeError(f"{src}: {sys.getdefaultencoding()} :{err}")


def expand_data(data_to_expand):
    """
    From data structure like setup.py data_files (see http://)
     [(directory/where/to/copy/the/file, [path/to/file/in/archive/file1, ...]), ...]
    but instead of the original struct this one accept to specify a directory in elements to copy.

    This function will generate one entry for all *content* of the directory and subdirectory
    recursively, to in fine copy the tree in archive in dest on the host

    the first level of directory itself is not include (which allow to rename it)
    :param data_to_expand:
    :type  data_to_expand: list of tuple
    :return: list of tuple
    """
    def remove_prefix(prefix, path):
        prefix = os.path.normpath(prefix)
        path = os.path.normpath(path)
        to_remove = len([i for i in prefix.split(os.path.sep) if i])
        truncated = [i for i in path.split(os.path.sep) if i][to_remove:]
        truncated = os.path.sep.join(truncated)
        return truncated

    data_struct = []

    for base_dest_dir, src in data_to_expand:
        base_dest_dir = os.path.normpath(base_dest_dir)
        for one_src in src:
            if os.path.isdir(one_src):
                for path, _, files in os.walk(one_src):
                    if not files:
                        continue
                    path_2_create = remove_prefix(one_src, path)
                    data_struct.append((os.path.join(base_dest_dir, path_2_create), [os.path.join(path, f) for f in files]))
            if os.path.isfile(one_src):
                data_struct.append((base_dest_dir, [one_src]))
    return data_struct


def read_md(f):
    import codecs
    with codecs.open(f, 'r', encoding='utf8') as f:
        text = f.read()
    return text


###################################################
# the configuration of the installer start bellow #
###################################################

setup(name='macsyfinder',
      version=mf_version,
      description="MacSyFinder: Detection of macromolecular systems in protein datasets using systems modelling and similarity search",
      long_description=read_md('README.md'),
      long_description_content_type="text/markdown",
      author="Sophie Abby, Bertrand Néron",
      author_email="sophie.abby@univ-grenoble-alpes.fr, bneron@pasteur.fr",
      url="https://github.com/gem-pasteur/macsyfinder/",
      download_url='https://github.com/gem-pasteur/macsyfinder/',
      license="GPLv3",
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
          ],
      python_requires='>=3.7',
      install_requires=[i for i in [l.strip() for l in open("requirements.txt").read().split('\n')] if i],
      extras_require={'dev': open("requirements_dev.txt").read().split()},
      test_suite='tests.run_tests.discover',
      zip_safe=False,
      packages=[p for p in find_packages() if p != 'tests'],
      # file where some variables must be fixed by install
      fix_prefix=['macsypy/__init__.py'],
      entry_points={
          'console_scripts': [
              'macsyfinder=macsypy.scripts.macsyfinder:main',
              'macsydata=macsypy.scripts.macsydata:main',
              'macsyprofile=macsypy.scripts.macsyprofile:main',
          ]
      },
      # (dataprefix +'where to put the data in the install, [where to find the data in the tar ball]
      data_files=expand_data([('share/macsyfinder/data/', ['data']),
                              ('share/doc/macsyfinder/html', ['doc/build/html']),
                              ('share/doc/macsyfinder/pdf', ['doc/build/latex/macsyfinder.pdf']),
                              ('etc/macsyfinder', ['etc/macsyfinder.conf'])
                              ]),

      cmdclass={'install_lib': install_lib},
      distclass=UsageDistribution
      )
