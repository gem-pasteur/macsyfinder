import sys
import os.path
import unittest


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