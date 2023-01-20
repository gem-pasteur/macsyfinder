.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2023  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _registries:


**********
registries
**********

The registry manage the different location where `macsyfinder` can find models definitions and their associated profiles.

.. _registries_api:

registries API reference
========================

ModelRegistry
=============
.. autoclass:: macsypy.registries.ModelRegistry
   :members:
   :private-members:
   :special-members:


ModelLocation
=============
.. autoclass:: macsypy.registries.ModelLocation
   :members:
   :private-members:
   :special-members:


MetaDefLoc
==========
.. autoclass:: macsypy.registries.MetaDefLoc
   :members:
   :private-members:
   :special-members:


DefinitionLocation
==================
.. autoclass:: macsypy.registries.DefinitionLocation
   :members:
   :private-members:
   :special-members:


split_def_name
===================
.. autofunction:: macsypy.registries.split_def_name


join_def_path
===================
.. autofunction:: macsypy.registries.join_def_path


scan_models_dir
===============
.. autofunction:: macsypy.registries.scan_models_dir