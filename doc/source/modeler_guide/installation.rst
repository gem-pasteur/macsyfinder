.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2024 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.


.. _modeler_installation:

************
Installation
************

MacSyFinder works with models for macromolecular systems that are not shipped with it,
you have to install them separately. See the :ref:`macsydata section <modeler_macsydata>` below.
We also provide container so you can use macsyfinder directly.


==================================
MacSyFinder Installation procedure
==================================

To develop new models and share them, MacSyFinder requires *git* and the python librarie *GitPython*

Below the procedure to install *MacSyFinder* in *modeler* mode in a virtualenv

.. code-block:: bash

    python3 -m venv macsyfinder
    cd macsyfinder
    source bin/activate
    python3 -m pip install macsyfinder[model]


*GitPython* dependency will be automatically retrieved and installed when using `pip` for installation (see below).

.. warning::

    But you have to install *git* by yourself (using your preferred package manager)


From Conda/Mamba
================

From version 2.0 and above, MacSyFinder is packaged for Conda/Mamba
The Conda/Mamba package include modeler dependencies


From container
==============

From version 2.0 and above, a docker image is available. This image allow you to develop models.


.. _modeler_macsydata:

====================================
Models installation with `macsydata`
====================================

Once MacSyFinder is installed you have access to an utility program to manage the models: `macsydata`

This script allows to search, download, install and get information from MacSyFinder models stored on
github (https://github.com/macsy-models) or locally installed. The general syntax for `macsydata` is::

    macsydata <general options> <subcommand> <sub command options> <arguments>


To list all models available on *macsy-models*::

    macsydata available

To search for models on *macsy-models*::

    macsydata search TXSS

you can also search in models description::

    macsydata search -S secretion

To install a model package::

    macsydata install <model name>

To install a model when you have not the right to install it system-wide

To install it in your home (*./macsyfinder/data*)::

    macsydata install --user <model name>

To install it in any directory::

    macsydata install --target <model dir> <model_name>

To know how to cite a model package::

    macsydata cite <model name>

To show the model definition::

    macsydata definition <package or subpackage> model1 [model2, ...]

for instance to show model definitions T6SSii and T6SSiii in TXSS+/bacterial subpackage::

    macsydata definition TXSS+/bacterial T6SSii T6SSiii

To show all models definitions in TXSS+/bacterial subpackage::

    macsydata definition TXSS+/bacterial

To create a git repository with a skeleton for your own model package::

    macsydata init --pack-name <MY_PACK_NAME> --maintainer <"mantainer name"> --email <maintainer email> --authors <"author1, author2, ..">

above macsydata with required options. Below I add optioanl but recommended options. ::

    macsydata init --pack-name <MY_PACK_NAME> --maintainer <"mantainer name"> --email <maintainer email> --authors <"author1, author2, .."> \
    --license cc-by-nc-sa --holders <"the copyright holders"> --desc <"one line package description">

To list all `macsydata` subcommands::

    macsydata --help

To list all available options for a subcommand::

    macsydata <subcommand> --help

For models not stored in *macsy-models* the commands *available*, *search*,
*installation* from remote or *upgrade* from remote are **NOT** available.

For models **NOT** stored in *macsy-models*, you have to manage them semi-manually.
Download the archive (do not unarchive it), then use *macsydata* to install the archive.
