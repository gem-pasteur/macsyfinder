.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  

.. _installation:


************
Installation
************

MacSyFinder works with models for macromolecular systems that are not shipped with it, 
you have to install them separately. See the :ref:`macsydata section <macsydata>` below.


.. _dependencies:

========================
MacSyFinder dependencies
========================
**Python version >=3.7** is required to run MacSyFinder: https://docs.python.org/3.7/index.html

MacSyFinder has one program dependency:

 - the *Hmmer* program, version 3.1 or greater (http://hmmer.org/).

The *hmmsearch* program should be installed (*e.g.*, in the PATH) in order to use MacSyFinder.
Otherwise, the paths to this executable must be specified in the command-line:
see the :ref:`command-line options <hmmer-options>`.
 
 
MacSyFinder also relies on four Python library dependencies:

 - colorlog
 - pyyaml
 - packaging
 - networkx

These dependencies will be automatically retrieved and installed when using `pip` for installation (see below). 
 

==================================
MacSyFinder Installation procedure
==================================

It is recommended to use `pip` to install the MacSyFinder package.

Archive overview
================

* **doc** => the documentation in html and pdf
* **etc** => a template of macsyfinder configuration file
* **test** => all what is needed for unitary tests
* **macsypy** => the macsyfinder python library
* **setup.py** => the installation script
* **requirements.txt** => the python dependencies
* **requirements_dev.txt** => the python dependencies for developers
* **COPYING** => the licensing
* **COPYRIGHT** => the copyright
* **README.md** => very brief macsyfinder overview
* **CONTRIBUTORS** => list of people who contributed to the code


Installation steps:
=======================

Make sure every required dependency/software is present.
--------------------------------------------------------

By default MacSyFinder will try to use `hmmsearch` in your PATH. If `hmmsearch` is not in the PATH,
you have to set the absolute path to `hmmsearch` in a :ref:`configuration file <config-definition-label>` 
or in the :ref:`command-line <hmmer-options>` upon execution.
If the tools are not in the path, some test will be skipped and a warning will be raised.


Perform the installation.
-------------------------

::

    pip install --no-binary macsyfinder macsyfinder


If you do not have the privileges to perform a system-wide installation,
you can either install it in your home directory or
use a `virtual environment <https://virtualenv.pypa.io/en/stable/>`_.

installation in your home directory
"""""""""""""""""""""""""""""""""""

::

    pip install --user --no-binary macsyfinder macsyfinder


installation in a virtualenv
""""""""""""""""""""""""""""

::

    python3.7 -m venv macsyfinder
    cd macsyfinder
    source bin/activate
    pip install --no-binary macsyfinder macsyfinder

To exit the virtualenv just execute the `deactivate` command.
To run `macsyfinder`, you need to activate the virtualenv: ::

    source macsyfinder/bin/activate

Then run `macsyfinder` or `macsydata`.

  
.. note::
  Super-user privileges (*i.e.*, ``sudo``) are necessary if you want to install the program in the general file architecture.
  
  
.. note::
  If you do not have the privileges, or if you do not want to install MacSyFinder in the Python libraries of your system, 
  you can install MacSyFinder in a virtual environment (http://www.virtualenv.org/).

.. warning::
  When installing a new version of MacSyFinder, do not forget to uninstall the previous version installed ! 


Uninstalling MacSyFinder
========================

To uninstall MacSyFinder (the last version installed), run::

  (sudo) pip uninstall macsyfinder

If you install it in a virtualenv, just delete the virtual environment.
For instance if you create a virtualenv name macsyfinder::

    python3.7 -m venv macsyfinder

To delete it, remove the directory::

    rm -R macsyfinder


.. _macsydata:

====================================
Models installation with `macsydata`
====================================

Once MacSyFinder is installed you have access to an utility program to manage the models: `macsydata`

This script allows to search, download, install and get information from MacSyFinder models stored on github or locally
installed. The general syntax for `macsydata` is::

    macsydata <general options> <subcommand> <sub command options> <arguments>


To list all models available::

    macsydata available

To search for models::

    macsydata search TXSS

you can also search in models description::

    macsydata search -S secretion

To install a model package::

    macsydata install <model name>

To install a model when you have not the right to install it system-wide::

    macsydata install --user <model name>

To know how to cite a model package::

    macsydata cite <model name>

To list all `macsydata` subcommands::

    macsydata --help

To list all available options for a subcommand::

    macsydata <subcommand> --help
