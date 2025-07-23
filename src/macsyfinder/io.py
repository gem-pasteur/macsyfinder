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

from macsylib import __version__ as msl_vers, __commit__ as msl_commit
from macsyfinder import __version__ as msf_vers, __commit__ as msf_commit

def outfile_header(models_fam_name: str,
                   models_version: str,
                   skipped_replicons: list[str] | None = None) -> str:
    """

    :return: The first lines of each result file
    """
    header = f"""# macsyfinder {msf_vers} {msf_commit} (using MacSyLib {msl_vers} {msl_commit})
# models : {models_fam_name}-{models_version}
# {' '.join(sys.argv)}"""
    if skipped_replicons:
        header += "\n#"
        for rep_name in skipped_replicons:
            header += f"\n# WARNING: The replicon '{rep_name}' has been SKIPPED. Cannot be solved before timeout."
        header += "\n#"
    return header
