.. _HMMReport:

*************
HMMReport API
*************

A *"HMMReport"* object represents the results of a Hmmer program search on a dataset with a hidden Markov model protein profile (see :ref:`this section <report-implementation>`).
This object has methods to extract and filter Hmmer raw outputs (see :ref:`generated output files <hmmer-outputs-label>`), and then build Hits relevant for system detection. 
For matches selected with the filtering parameters, *"Hit"* objects (:class:`macsypy.HMMReport.Hit`) are built. 



HMMReport API reference
=======================

.. automodule:: macsypy.report
   :members: HMMReport
   :private-members:
   :special-members:

GeneralHMMReport API reference
==============================

.. automodule:: macsypy.report
   :members: GeneralHMMReport
   :private-members:
   :special-members:

OrderedHMMReport
================

.. automodule:: macsypy.report
   :members: OrderedHMMReport
   :private-members:
   :special-members:
   
GembaseHMMReport 
================

.. automodule:: macsypy.report
   :members: GembaseHMMReport
   :private-members:
   :special-members:

.. _hit-label:

Hit 
===

.. automodule:: macsypy.report
   :members: Hit
   :private-members:
   :special-members:
