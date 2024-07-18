.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2023  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _profile:

*******
profile
*******

The :ref:`Profile object <profile-implementation>` is used for the search of the gene with Hmmer.
A *"Profile"* must match a HMM protein profile file, which name is based on the profile name.
For instance, the *gspG* gene has the corresponding "gspG.hmm" profile file provided at a dedicated location.

.. _profile_api:

profile API reference
=====================

ProfileFactory
==============
 .. autoclass:: macsypy.profile.ProfileFactory
   :members:
   :private-members:
   :special-members:


Profile
=======

.. autoclass:: macsypy.profile.Profile
   :members:
   :private-members:
   :special-members:
