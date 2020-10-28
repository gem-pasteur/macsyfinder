.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.
    
.. _HMMReport:

******
report
******

A *"HMMReport"* object represents the results of a Hmmer program search on a dataset with a hidden Markov model protein profile (see :ref:`this section <report-implementation>`).
This object has methods to extract and filter Hmmer raw outputs (see :ref:`generated output files <hmmer-outputs-label>`), and then build Hits relevant for system detection. 
For matches selected with the filtering parameters, *"Hit"* objects (:class:`macsypy.HMMReport.Hit`) are built. 



HMMReport API reference
=======================

.. autoclass:: macsypy.report.HMMReport
   :members:
   :private-members:
   :special-members:

GeneralHMMReport API reference
==============================

.. autoclass:: macsypy.report.GeneralHMMReport
   :members:
   :private-members:
   :special-members:

OrderedHMMReport
================

.. autoclass:: macsypy.report.OrderedHMMReport
   :members:
   :private-members:
   :special-members:
   
GembaseHMMReport 
================

.. autoclass:: macsypy.report.GembaseHMMReport
   :members:
   :private-members:
   :special-members:


