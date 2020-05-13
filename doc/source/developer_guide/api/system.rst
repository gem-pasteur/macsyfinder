.. MacSyFinder - Detection of macromolecular systems in protein datasets
   using systems modelling and similarity search.
   Authors: Sophie Abby, Bertrand Néron
   Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
   See the COPYRIGHT file for details
   MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
   See the COPYING file for details.
    
.. _system:

******
System
******

It represents an occurrence of a model in a replicon.


System
=======
.. autoclass:: macsypy.system.System
   :members:
   :private-members:
   :special-members:

SystemSerializer
================

Use to serilize the systems found. :class:`macsypy.system.SystemSerializer` is an abstract class.
There is 3 concretes classes which implements SystemSerializer:

    * :class:`macsypy.system.TxtSystemSerializer`
    * :class:`macsypy.system.TsvSystemSerializer`
    * :class:`macsypy.system.JsonSystemSerializer`

.. automodule:: macsypy.system
   :members: SystemSerializer, TxtSystemSerializer, TsvSystemSerializer, JsonSystemSerializer
   :private-members:
   :special-members: