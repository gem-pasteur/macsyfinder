.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _package:


*******
package
*******

Allow to handles model package either on localhost or from a remote location.
the model packages can be stored in github organization to be downloaded and installed locally.
The classes below are used by `macsydata`, which is the entry point to manipulate models package.

.. _model_index_api:

ModelIndex
==========
.. automodule:: macsypy.package
   :members: AbstractModelIndex, LocalModelIndex, RemoteModelIndex
   :private-members:
   :special-members:


.. _package_api:

Package
=======
.. autoclass:: macsypy.package.Package
   :members:
   :private-members:
   :special-members:

