.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _publish_package:

*************************
Publishing/sharing models
*************************


 .. _model_package:

Writing my own macsy-model package
==================================

The whole package structure is described :ref:`above <package_structure>` and requires five different types of files described below to be complete:

* a metadata file
* a README.md file
* a LICENCE file
* macsy-models definition(s)
* HMM profiles


metadata file
-------------

This file contains some meta information about the package itself.
It is in `YAML <https://en.wikipedia.org/wiki/YAML>`_ format and must have the following structure:

.. code-block:: yaml

    ---
    maintainer:
      name: The name of the person who maintains/to contact for further information. (required)
      email: The email of the maintainer (required)
    short_desc: A one line description of the package (can e.g. be used for *macsydata* searches) (required)
    vers: The package version (required)
    cite: The publication(s) to cite by the user when the package is used (optional, used by `macsydata cite`)
    doc: Where to find extended documentation (optional)
    licence: The licence under the package is released (optional but highly recommended)
    copyright: The copyright of the package (optional)

For example:

.. code-block:: yaml

    ---
    maintainer:
       name: first name last name
       email: login@my_domain.com
    short_desc: Models for 15 types of secretion systems or bacterial appendages (T1SS, T2SS, T3SS, T4P, pT4SSt, pT4SSi, T5aSS, T5bSS, T5bSS, T6SSi, T6SSii, T6SSiii, Flagellum, Tad, T9SS).
    vers: 0.0a1
    cite:
       - |
         Abby Sophie S., Cury Jean, Guglielmini Julien, Néron Bertrand, Touchon Marie, Rocha Eduardo P. C. (2016).
         Identification of protein secretion systems in bacterial genomes.
         In Scientific Reports, 6, pp. 23080.
         http://dx.doi.org/10.1038/srep23080
    doc: https://github.com/macsy-models/TXSS
    licence: CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
    copyright: 2014-2020, Institut Pasteur, CNRS

.. warning::
    This `metadata.yml` file is **mandatory**. Without this file your archive/repository will not be considered as a *macsy-model package*.

.. note::

    * *-* specify an item of yaml list
    * *|* is used to specify a single item but over multiple lines.



Share your models
=================

If you want to share your models you can create a :ref:`macsy-model package <model_package>` in your github repository.

1. check the validity of your package with the ``macsydata check`` command.
2. create a tag, and submit a pull request to https://github.com/macsy-models organization.
3. when your pull request (PR) is accepted, the model package becomes automatically available to the community through the *macsydata* tool.

If you don't want to submit a PR you can provide the tag release tarball (tar.gz) as is to your collaborators.
This archive will also be usable with the `macsydata` tool.
