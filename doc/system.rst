.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014  Institut Pasteur, Paris.                           
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _system:

**********
System API
**********

It represents a macromolecular system to detect. See :ref:`system-implementation` for an overview of its implementation.

.. note:: 
    
    a complete description of macromolecular systems modelling is available in the section :ref:`system_definition`
    

SystemBank API reference
========================
 .. automodule:: macsypy.system
   :members: SystemBank
   :private-members:
   :special-members:

 
.. note::

   Don't instanciate your own SystemBank use the system_bank at the top level of the module. ::
     
     from macsypy.system import system_bank
 
System API reference
====================

.. automodule:: macsypy.system
   :members: System
   :private-members:
   :special-members:
