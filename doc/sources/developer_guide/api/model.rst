.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014  Institut Pasteur, Paris.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _model:


*********
model API
*********




SystemBank API reference
========================
 .. automodule:: macsypy.model
   :members: ModelBank
   :private-members:
   :special-members:


.. note::

   Don't instanciate your own SystemBank use the system_bank at the top level of the module. ::

     from macsypy.system import system_bank

System API reference
====================

.. automodule:: macsypy.model
   :members: Model
   :private-members:
   :special-members: