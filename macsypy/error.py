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


class MacsypyError(Exception):
    """
    The base class for MacSyFinder specific exceptions.
    """
    pass


class MacsydataError(MacsypyError):
    """
    Raised when error is encounter during model package handling
    """
    pass


class MacsyDataLimitError(MacsydataError):
    """
    Raised when the maximum number of github api call is reached
    """
    pass


class OptionError(MacsypyError):
    """
    Raised when command line option is not set properly
    """
    pass


class ModelInconsistencyError(MacsypyError):
    """
    Raised when a definition model is not consistent.
    """
    pass


class SystemDetectionError(MacsypyError):
    """
    Raised when the detection of systems from Hits encountered a problem.
    """
    pass
