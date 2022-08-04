.. MacSyFinder - Detection of macromolecular systems in protein datasets
   using systems modelling and similarity search.
   Authors: Sophie Abby, Bertrand Néron
   Copyright © 2014-2022  Institut Pasteur (Paris), and CNRS.
   See the COPYRIGHT file for details
   MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
   See the COPYING file for details.
    
.. _serialization:

*************
serialization
*************

This module is a technical module where we can find the different way
to serialize the results:

   * the Systems found
   * The best solutions (best combination of systems)
   * The rejected candidates


.. _serialization_api:

SystemSerializer
================
.. autoclass:: macsypy.serialization.SystemSerializer
   :members:
   :private-members:
   :special-members:


TsvSystemSerializer
===================
.. autoclass:: macsypy.serialization.TsvSystemSerializer
   :members:
   :private-members:
   :special-members:


TsvSolutionSerializer
=====================
.. autoclass:: macsypy.serialization.TsvSolutionSerializer
   :members:
   :private-members:
   :special-members:


TsvLikelySystemSerializer
=========================
.. autoclass:: macsypy.serialization.TsvLikelySystemSerializer
   :members:
   :private-members:
   :special-members:


TsvRejectedCandidatesSerializer
===============================
.. autoclass:: macsypy.serialization.TsvRejectedCandidatesSerializer
   :members:
   :private-members:
   :special-members:


TsvSpecialHitSerializer
=======================
.. autoclass:: macsypy.serialization.TsvSpecialHitSerializer
   :members:
   :private-members:
   :special-members:


TxtSystemSerializer
===================
.. autoclass:: macsypy.serialization.TxtSystemSerializer
   :members:
   :private-members:
   :special-members:


TxtLikelySystemSerializer
=========================
.. autoclass:: macsypy.serialization.TxtLikelySystemSerializer
   :members:
   :private-members:
   :special-members:


TxtUnikelySystemSerializer
==========================
.. autoclass:: macsypy.serialization.TxtUnikelySystemSerializer
   :members:
   :private-members:
   :special-members:





