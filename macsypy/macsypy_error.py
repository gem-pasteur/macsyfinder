# -*- coding: utf-8 -*-

###################################
# Created on Jul 15, 2013
#
# @author: bneron
# @contact: user_email
# @organization: Institut Pasteur
# @license: license
###################################


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
