import os
import sys
import unittest
import platform
from StringIO import StringIO
from contextlib import contextmanager
import hashlib
from functools import partial


class MacsyTest(unittest.TestCase):

    _data_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "data"))

    def setsid(self):
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


    def assertFileEqual(self, f1, f2, msg=None):
        self.maxDiff = None
        with open(f1) as fh1, open(f2) as fh2:
            self.assertMultiLineEqual(fh1.read(), fh2.read(), msg=msg)

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


def md5sum(filename):
    """Compute file md5 checksum."""
    with open(filename, mode='rb') as f:
        d = hashlib.md5()
        for buf in iter(partial(f.read, 128), b''):
            d.update(buf)
    return d.hexdigest()
