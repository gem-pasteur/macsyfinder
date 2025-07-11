##########################################################################
#  MacSyFinder - Detection of macromolecular systems in protein dataset  #
#                using systems modelling and similarity search.          #
#  Authors: Sophie Abby, Bertrand Neron                                  #
#  Copyright (c) 2014-2025  Institut Pasteur (Paris) and CNRS.           #
#  See the COPYRIGHT file for details                                    #
#                                                                        #
#  This file is part of MacSyFinder package.                             #
#                                                                        #
#  MacSyFinder is free software: you can redistribute it and/or modify   #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  MacSyFinder is distributed in the hope that it will be useful,        #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details .                         #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with MacSyFinder (COPYING).                                     #
#  If not, see <https://www.gnu.org/licenses/>.                          #
##########################################################################

import sys
from macsyfinder import __version__ as msf_ver, __commit__ as msf_commit
from macsylib import __version__ as msl_ver, __commit__ as msl_commit
from macsyfinder.io import outfile_header
from . import MacsyTest


class IoTest(MacsyTest):


    def test_oufile_header(self,):
        model_fam = 'models'
        model_ver = '0.0'
        expected = f"""# macsyfinder {msf_ver} {msf_commit} (using MacSyLib {msl_ver} {msl_commit})
# models : models-0.0
# {' '.join(sys.argv)}"""
        self.assertEqual(expected, outfile_header(model_fam, model_ver))

        expected = f"""# macsyfinder {msf_ver} {msf_commit} (using MacSyLib {msl_ver} {msl_commit})
# models : models-0.0
# {' '.join(sys.argv)}
#
# WARNING: The replicon 'rep_1' has been SKIPPED. Cannot be solved before timeout.
# WARNING: The replicon 'rep_2' has been SKIPPED. Cannot be solved before timeout.
#"""

        self.assertEqual(expected,
                         outfile_header(model_fam, model_ver,
                                               ['rep_1', 'rep_2']))
