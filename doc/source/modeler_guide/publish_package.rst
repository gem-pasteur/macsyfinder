.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2021 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _publish_package:

*************************
Publishing/sharing models
*************************


 .. _writing_model_package:

Writing my own macsy-model package
==================================

The whole package structure is described :ref:`above <package_structure>` and requires five different
types of files described below to be complete:

* a metadata.yml file (mandatory)
* a README.md file (mandatory)
* a LICENCE file (optional but **HIGHLY** recommended)
* a model_conf.xml file (optional)
* macsy-models definition(s) (mandatory)
* HMM profiles (mandatory)

Writing metadata
----------------

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
    license: The license under the package is released (optional but highly recommended)
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
    license: CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
    copyright: 2014-2021, Institut Pasteur, CNRS

.. warning::
    This `metadata.yml` file is **mandatory**. Without this file your archive/repository will not be considered as a *macsy-model package*.

.. note::

    * *-* specify an item of yaml list
    * *|* is used to specify a single item but over multiple lines.



Share your models
=================

If you want to share your models you can create a :ref:`macsy-model package <model_package>` in your github repository.

1. check the validity of your package with the ``macsydata check`` command.
   macsydata check will report you:

   * every thing is clear: macsydata display you the next step to perform to publish the package
   * warning : it mean that the package could be improved.
             it better to fix it, but you go to *step 3*
   * error: the package is not ready to be publish as is. fix the errors before to go to *step 3*

2. create a tag, and submit a pull request to https://github.com/macsy-models organization.
   This step is **very important**, without tag, there is no package.
   macsydata check only tagged packages.
   It's also of the model provider to setup a tag with the same name as the version in the metadata file.
   It is recommended to follow a versioning scheme describe here https://www.python.org/dev/peps/pep-0440/#public-version-identifiers
3. when your pull request (PR) is accepted, the model package becomes automatically available to the community through the *macsydata* tool.

If you don't want to submit a PR you can provide the tag release tarball (tar.gz) as is to your collaborators.
This archive will also be usable with the `macsydata` tool.

.. note:: macsydata check
    check the syntax of the package. But it does not publish anything.
    It just warn you if something goes wrong with the package.
    every model provider should check it's own package before to publish it.
    The package publication is done by the `git push` and the `pull request`.

example of macsydata check outputs

your package is syntactically correct

.. code-block:: text

    macsydata check tests/data/models/test_model_package/
    Checking 'test_model_package' package structure
    Checking 'test_model_package' metadata_path
    Checking 'test_model_package' Model definitions
    Models Parsing
    Definitions are consistent
    Checking 'test_model_package' model configuration
    There is no model configuration for package test_model_package.
    If everyone were like you, I'd be out of business
    To push the models in organization:
            cd tests/data/models/test_model_package
    Transform the models into a git repository
            git init .
            git add .
            git commit -m 'initial commit'
    add a remote repository to host the models
    for instance if you want to add the models to 'macsy-models'
            git remote add origin https://github.com/macsy-models/
            git tag 1.0b2
            git push --tags


You received some warnings

.. code-block:: text

    macsydata check tests/data/models/Model_w_conf/
    Checking 'Model_w_conf' package structure
    Checking 'Model_w_conf' metadata_path
    Checking 'Model_w_conf' Model definitions
    Models Parsing
    Definitions are consistent
    Checking 'Model_w_conf' model configuration
    The package 'Model_w_conf' have not any LICENSE file. May be you have not right to use it.
    The package 'Model_w_conf' have not any README file.
    macsydata says: You're only giving me a partial QA payment?
    I'll take it this time, but I'm not happy.
    I'll be really happy, if you fix warnings above, before to publish these models.

You received some errors

.. code-block:: text

    macsydata check tests/data/models/TFF-SF/
    Checking 'TFF-SF' package structure
    The package 'TFF-SF' have no 'metadata.yml'.
    Please fix issues above, before publishing these models.
    ValueError
