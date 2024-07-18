.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2023  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _configuration:

*************
configuration
*************

Options to run MacSyFinder can be specified in a :ref:`Configuration file <config-definition-label>`.
The API described below handles all configuration options for MacSyFinder.

The :class:`macsypy.config.MacsyDefaults` hold the default values for `macsyfinder` whereas
the :class:`macsypy.config.Config` hold the values for a `macsyfinder` run.

.. _config_api:

configuration API reference
===========================

MacsyDefaults
=============
.. autoclass:: macsypy.config.MacsyDefaults
   :members:
   :private-members:
   :special-members:


Config
======
.. autoclass:: macsypy.config.Config
   :members:
   :private-members:
   :special-members:


NoneConfig
==========
.. autoclass:: macsypy.config.NoneConfig
   :members:
   :private-members:
   :special-members:
