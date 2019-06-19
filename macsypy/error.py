# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


class MacsypyError(Exception):
    """
    The base class for MacSyFinder specific exceptions.
    """
    pass

class OptionError(MacsypyError):
    """
    Raised when command line option is not set properly
    """
    pass

class ModelInconsistencyError(MacsypyError):
    """
    Raised when a secretion model is not consistent.
    """
    pass


class SystemDetectionError(MacsypyError):
    """
    Raised when the detection of systems from Hits encountered a problem.
    """
    pass
