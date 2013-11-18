.. _HMMReport:

*************
HMMReport API
*************

A *"HMMReport"* object represents the results of a Hmmer program search on a dataset with a hidden Markov model protein profile (see :ref:`this section <report-implementation>`).
This object has methods to extract and filter Hmmer raw outputs (see :ref:`generated output files <hmmer-outputs-label>`), and then build Hits relevant for system detection. 
For matches selected with the filtering parameters, *"Hit"* objects (:class:`txsscanlib.HMMReport.Hit`) are built. 



HMMReport API reference
=======================

.. automodule:: txsscanlib.report
   :members: HMMReport
   :private-members:
   :special-members:

GeneralHMMReport API reference
==============================

.. automodule:: txsscanlib.report
   :members: GeneralHMMReport
   :private-members:
   :special-members:

OrderedHMMReport
================

.. automodule:: txsscanlib.report
   :members: OrderedHMMReport
   :private-members:
   :special-members:
   
GembaseHMMReport 
================

.. automodule:: txsscanlib.report
   :members: GembaseHMMReport
   :private-members:
   :special-members:

.. _hit-label:

Hit 
===

.. automodule:: txsscanlib.report
   :members: Hit
   :private-members:
   :special-members:
