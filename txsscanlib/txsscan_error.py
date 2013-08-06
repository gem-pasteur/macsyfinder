# -*- coding: utf-8 -*-

###################################
# Created on Jul 15, 2013
#
# @author: bneron
# @contact: user_email
# @organization: Institut Pasteur
# @license: license
###################################


class TxsscanError(Exception):
    """
    The base class for txsscan specific exceptions.
    """
    pass


class SystemInconsistencyError(TxsscanError):
    """
    Raised when a secretion system is not consistent.
    """
    pass

class SystemDetectionError(TxsscanError):
    """
    Raised when the detection of systems from Hits encountered a problem.
    """
    pass
