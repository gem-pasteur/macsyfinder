#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2025  Institut Pasteur (Paris) and CNRS.           #
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
"""
MacSyFinder package contains mainly variable used in library as __version__, __citation__
and helper functions
"""

import os
import sys
import subprocess


__version__ = '2.1.6'


__citation__ = """Néron, Bertrand; Denise, Rémi; Coluzzi, Charles; Touchon, Marie; Rocha, Eduardo P.C.; Abby, Sophie S.
MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes.
Peer Community Journal, Volume 3 (2023), article no. e28. doi : 10.24072/pcjournal.250.
https://peercommunityjournal.org/articles/10.24072/pcjournal.250/"""


def get_git_revision_short_hash() -> str:
    """
    :return: the git commit number (short version) or empty string if this not a git repository
    :rtype: str
    """
    try:
        short_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'],
                                             cwd=os.path.dirname(os.path.abspath(__file__)),
                                             stderr=subprocess.DEVNULL)

        short_hash = str(short_hash, "utf-8").strip()
    except Exception:
        short_hash = ''
    return short_hash

#  do not display the commit for the MSF tagged versions
__commit__ = f'{get_git_revision_short_hash()}' if 'dev' in __version__ else ''


def get_version_message() -> str:
    """
    :return: the long description of the macsyfinder version
    """
    # I do not understand why
    # BUT if i do the import at the top level
    # pip cannot install macsyfinder with an obscured error
    # ValueError: malformed node or string on line 30: <ast.JoinedStr object at 0x7fe8cc127c10>
    # ...
    # AttributeError: macsyfinder has no attribute __version__
    from macsylib import __version__ as msl_version, __commit__ as msl_commit  # pylint: disable=import-outside-toplevel
    from macsylib import solution  # pylint: disable=import-outside-toplevel
    import pandas as pd  # pylint: disable=import-outside-toplevel

    version = __version__
    commit = __commit__
    py_vers = sys.version.replace('\n', ' ')
    vers_msg = f"""MacSyFinder {version} {commit}
using:
- Python {py_vers}
- MacSyLib {msl_version} {msl_commit}
- NetworkX {solution.nx.__version__}
- Pandas {pd.__version__}

MacSyFinder is distributed under the terms of the GNU General Public License (GPLv3).
See the COPYING file for details.

If you use this software please cite:
{__citation__}

and don't forget to cite the used models:
msf_data cite <model>
"""
    return vers_msg
