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
wrapper to macsylib.macsyprofile
"""
import warnings
import textwrap

from macsylib.scripts.macsyprofile import main as msl_main

def _cmde_line_header():
    return textwrap.dedent(r'''

         *            *                 *              *
    *           *                 *   *   *  *    **
       *   *         __       *         **    __ _ _     *
      _ __ ___  ___ / _|     _ __  _ __ ___  / _(_) | ___
     | '_ ` _ \/ __| |_     | '_ \| '__/ _ \| |_| | |/ _ \
     | | | | | \__ \  _|    | |_) | | | (_) |  _| | |  __/
     |_| |_| |_|___/_|______| .__/|_|  \___/|_| |_|_|\___|
          *    *    |_______|_|     *               **
     *      *   * *     *   **         *   *  *
      *      *         *        *    *              *
    *                           *  *           *        *

    msf_profile - Profile Helper Tool
    ''')

def deprecated_main(args: list[str] = None) -> None:
    warnings.warn(r"""
       __        __               _
       \ \      / /_ _ _ __ _ __ (_)_ __   __ _
        \ \ /\ / / _` | '__| '_ \| | '_ \ / _` |
         \ V  V / (_| | |  | | | | | | | | (_| |
          \_/\_/ \__,_|_|  |_| |_|_|_| |_|\__, |
                                          |___/

'macsyprofile' has been renamed. Use 'msf_profile' instead.
'macsyprofile' will be removed in future versions.""", DeprecationWarning, stacklevel=2)
    msl_main(args=args, package_name='macsyfinder', tool_name='msf_profile')


def main(args: list[str] = None) -> None:
    ############## WORKAROUND #####################
    # due to a bug in macsylib 1.0.4
    # I need to monkey patch MacsyDefaults
    import macsylib.scripts.macsyprofile
    # I override the default values of MacsyDefaults.__init__
    macsylib.scripts.macsyprofile.MacsyDefaults.__init__.__defaults__ = (None, 'macsyfinder')
    ###########  END of WORKAROUND  ###############
    msl_main(args=args,
             header=_cmde_line_header(),
             package_name='macsyfinder', tool_name='msf_profile')


if __name__ == "__main__":
    main()
