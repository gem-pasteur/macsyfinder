.. _system:

**********
System API
**********

It represents a macromolecular system to detect. See :ref:`system-implementation` for an overview of its implementation.

.. note:: 
    
    a complete description of the secretion system modelling is available in the section :ref:`system_definition`
    

SystemBank API reference
========================
 .. automodule:: txsscanlib.system
   :members: SystemBank
   :private-members:
   :special-members:

 
.. note::

   Don't instanciate your own SystemBank use the system_bank at the top level of the module. ::
     
     from txsscanlib.system import system_bank
 
System API reference
====================

.. automodule:: txsscanlib.system
   :members: System
   :private-members:
   :special-members:
