# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur (Paris) and CNRS.                                   #
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


class SystemInconsistencyError(MacsypyError):
    """
    Raised when a secretion system is not consistent.
    """
    pass


class SystemDetectionError(MacsypyError):
    """
    Raised when the detection of systems from Hits encountered a problem.
    """
    pass
