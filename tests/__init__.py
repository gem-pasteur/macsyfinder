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
import sys
import shutil
import unittest
import platform
from io import StringIO
from contextlib import contextmanager
import hashlib
from functools import partial
import tempfile
import uuid
import colorlog
import json

import macsypy
import macsypy.config


def path_to_modulename(p):
    """
    Example: foo/bar.py become bar
    """
    filename = os.path.basename(p)
    modulename = os.path.splitext(filename)[0]
    return modulename


class MacsyTest(unittest.TestCase):

    _tests_dir = os.path.normpath(os.path.dirname(__file__))
    _data_dir = os.path.join(_tests_dir, "data")

    def __init__(self, *args, **kwargs):
        macsypy.__MACSY_DATA__ = self._tests_dir
        macsypy.config.__MACSY_DATA__ = self._tests_dir
        super().__init__(*args, **kwargs)

    @staticmethod
    def setsid():
        platform = sys.platform
        if platform.startswith('linux'):
            setsid = 'setsid'
        elif platform.startswith('darwin'):
            setsid = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'utils', 'setsid'))
        else:
            setsid = ''
        return setsid

    @classmethod
    def find_data(cls, *args):
        data_path = os.path.join(cls._data_dir, *args)
        if os.path.exists(data_path):
            return data_path
        else:
            raise IOError("data '{}' does not exists".format(data_path))


    @contextmanager
    def catch_io(self, out=False, err=False):
        """
        Catch stderr and stdout of the code running within this block.
        """
        old_out = sys.stdout
        new_out = old_out
        old_err = sys.stderr
        new_err = old_err
        if out:
            new_out = StringIO()
        if err:
            new_err = StringIO()
        try:
            sys.stdout, sys.stderr = new_out, new_err
            yield sys.stdout, sys.stderr
        finally:
            sys.stdout, sys.stderr = old_out, old_err

    @staticmethod
    def fake_exit(*args, **kwargs):
        returncode = args[0]
        raise TypeError(returncode)

    @staticmethod
    def mute_call(call_ori):
        """
        hmmsearch or prodigal write lot of things on stderr or stdout
        which noise the unit test output
        So I replace the `call` function in module integron_finder
        by a wrapper which call the original function but add redirect stderr and stdout
        in dev_null
        :return: wrapper around call function
        :rtype: function
        """
        def wrapper(*args, **kwargs):
            with open(os.devnull, 'w') as f:
                kwargs['stderr'] = f
                kwargs['stdout'] = f
                res = call_ori(*args, **kwargs)
            return res
        return wrapper

    def assertFileEqual(self, f1, f2, comment=None, msg=None):
        self.maxDiff = None
        # the StringIO does not support context in python2.7
        # so we can use the following statement only in python3
        from itertools import zip_longest
        with open(f1) if isinstance(f1, str) else f1 as fh1, open(f2) if isinstance(f2, str) else f2 as fh2:
            for l1, l2 in zip_longest(fh1, fh2):
                if l1 and l2:
                    if comment and l1.startswith(comment) and l2.startswith(comment):
                        continue
                    self.assertEqual(l1, l2, msg)
                elif l1:  # and not l2
                    raise self.failureException(f"{fh1.name} is longer than {fh2.name}")
                elif l2:  # and not l1
                    raise self.failureException(f"{fh2.name} is longer than {fh1.name}")

    def assertTsvEqual(self, f1, f2, comment="#", msg=None):
        # the StringIO does not support context in python2.7
        # so we can use the following statement only in python3
        from itertools import zip_longest
        with open(f1) if isinstance(f1, str) else f1 as fh1, open(f2) if isinstance(f2, str) else f2 as fh2:
            header = True
            for i, grp in enumerate(zip_longest(fh1, fh2), 1):
                l1, l2 = grp
                if l1.startswith(comment) and l2.startswith(comment):
                    continue
                fields_1 = l1.split()
                fields_2 = l2.split()
                if not fields_1 and not fields_2:
                    # skip empty line
                    continue

                # the system_id may change from one run to another
                # So I have to remove them before to compare each row
                if header:
                    fields_nb = len(fields_1)
                    header = False
                if fields_nb == 20:  # all_systems.tsv
                    fields_1.pop(5)
                    fields_2.pop(5)
                elif fields_nb == 21:  # all_best_systems best_systems
                    fields_1.pop(6)
                    fields_2.pop(6)
                else:
                    raise RuntimeError(f"{fh1.name} {len(fields_1)}")

                if len(fields_1) == fields_nb - 1:
                    # remove used_in field if present
                    fields_1.pop(-1)
                    fields_2.pop(-1)
                self.assertListEqual(fields_1, fields_2, f"{fh1.name} differ from {fh2.name} at line {i}:\n{l1}{l2}")


    def assertSeqRecordEqual(self, s1, s2):
        for attr in ('id', 'name', 'seq'):
            s1_attr = getattr(s1, attr)
            s2_attr = getattr(s2, attr)
            self.assertEqual(s1_attr, s2_attr, msg="{} are different: {} != {}".format(attr, s1_attr, s2_attr))

        # there is a bug in some biopython version
        self.assertEqual(s1.description.rstrip('.'), s2.description.rstrip('.'))
        for s1_feat, s2_feat in zip(s1.features, s2.features):
            # location cannot be directly compared
            self.assertEqual(str(s1_feat.location), str(s2_feat.location))

            for attr in ('qualifiers', 'strand', 'type'):
                f1_attr = getattr(s1_feat, attr)
                f2_attr = getattr(s2_feat, attr)
                self.assertEqual(f1_attr, f2_attr, msg="{} are different: {} != {}".format(attr, f1_attr, f2_attr))

    def assertHmmEqual(self, hmm1, hmm2):
        with open(hmm1) as hmm1_file, open(hmm2) as hmm2_file:
            for hmm1_line, hmm2_line in zip(hmm1_file, hmm2_file):
                if hmm1_line.startswith('#') and hmm2_line.startswith('#'):
                    continue
                hmm1_fields = hmm1_line.split('#')[:-1]
                hmm2_fields = hmm2_line.split('#')[:-1]
                self.assertListEqual(hmm1_fields, hmm2_fields)


    def assertJsonEqual(self, json_file_1, json_file_2, max_diff=640):
        with open(json_file_1) as f1:
            j1 = json.load(f1)
        with open(json_file_2) as f2:
            j2 = json.load(f2)

        self.maxDiff = max_diff
        self.assertListEqual(j1, j2)


    @staticmethod
    def get_tmp_dir_name():
        return os.path.join(tempfile.gettempdir(), "macsyfinder_test_run")

    @staticmethod
    def get_uniq_tmp_dir_name():
        return os.path.join(tempfile.gettempdir(), "macsyfinder-{}".format(uuid.uuid4()))

    @staticmethod
    def rmtree(path):
        """
        Remove directory tree.

        :param path: the path to remove
        :type path: str
        """
        try:
            shutil.rmtree(path)
        except:
            pass

    @staticmethod
    def md5sum(file_=None, str_=None):
        """Compute md5 checksum.

        :param file_: the name of the file to compute the checksum for
        :type file_: str
        :param str_: the string to compute the checksum for
        :type str_: str
        """
        assert not (file_ and str_)

        d = hashlib.md5()

        if file_:
            with open(file_, mode='rb') as f:
                for buf in iter(partial(f.read, 128), b''):
                    d.update(buf)
        elif str_:
            assert isinstance(str_, str)
            d.update(str_)
        else:
            assert False
        return d.hexdigest()


    @contextmanager
    def catch_log(self, log_name='macsypy'):
        logger = colorlog.getLogger(log_name)
        handlers_ori = logger.handlers
        fake_handler = colorlog.StreamHandler(StringIO())
        try:
            logger.handlers = [fake_handler]
            yield LoggerWrapper(logger)
        finally:
            logger.handlers = handlers_ori


class LoggerWrapper(object):

    def __init__(self, logger):
        self.logger = logger

    def __getattr__(self, item):
        return getattr(self.logger, item)

    def get_value(self):
        return self.logger.handlers[0].stream.getvalue()


def which(name, flags=os.X_OK):
    """
    Search PATH for executable files with the given name.

    :param name: the name of the executable to search
    :type name: str
    :param flags: os mod the name must have, default is executable (os.X_OK).
    :type flags: os file mode R_OK|R_OK|W_OK|X_OK
    :return: the path of the executable
    :rtype: string or None
    """
    result = None
    path = os.environ.get('PATH', None)
    if path is None:
        return result
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if platform.system() == 'Windows':
            p += '.exe'
        if os.access(p, flags):
            result = p
            break
    return result
