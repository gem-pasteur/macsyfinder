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

MacSyFinder works with models, the models are not shipped with MacSyFinder.
So you have to install them separately.

========================
MacSyFinder dependencies
========================
MacSyFinder has two dependencies:
 - the program *Hmmer* version 3.1 (http://hmmer.janelia.org/).

*hmmsearch* programs should be installed (*e.g.*, in the PATH) in order to use MacSyFinder.
Otherwise, the paths to these executables must be specified in the command-line:
see the :ref:`command-line options <hmmer-options>`.
 
**Python version >=3.7** is required to run MacSyFinder: https://docs.python.org/3.7/index.html

==================================
MacSyFinder Installation procedure
==================================

It is recommended to use pip to install MacSyFinder package.

Archive overview
================

* **doc** => the documentation in html and pdf
* **etc** => a template of macsyfinder configuration file
* **test** => all needed for unit tests
* **macsypy** => the macsyfinder python library
* **setup.py** => the installation script
* **requirements.txt** => the python dependencies
* **requirements_dev.txt** => the python dependencies for developers
* **COPYING** => the licensing
* **COPYRIGHT** => the copyright
* **README.md** => very brief macsyfinder overview
* **CONTRIBUTORS** => list of people who contribute to the code


Installation steps:
=======================

Make sure every required dependence/software is present.
--------------------------------------------------------

By default MacSyFinder will try to use `hmmsearch` in your PATH, if `hmmsearch` is not in the PATH,
you have to set the absolute path to the `hmmsearch` in configuration.
If the tools are not in the path, some test will be skipped and a warning will be raised.


Perform the installation.
-------------------------

::

    pip install --no-binary macsyfinder macsyfinder


If you have not the privileges to perform a system wide installation,
you can use either install it in your home or
use a [virtual environment](https://virtualenv.pypa.io/en/stable/).

installation in your home
"""""""""""""""""""""""""

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

Then run macsyfinder/macsydata

  
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

If you install it in a virtualenv just delete the virtualenv


===================
Models installation
===================

Once MacSyFinder is installed you have access to an utility program to manage the models `macsydata`

this script allow to search, download, install get information from models store on github or from you locally
instaled models. The general syntax for `macsypy` is::

    macsydata <general options> <subcommand> <sub command options> <arguments>


for instance to list all models available::

    macsydata available

to search models::

    macsydata search TXSS

you can also search in models description::

    macsydata search -S secretion

To install a model package::

    macsydata install <model name>

To install a model when you have not the right to install in system wide::

    macsydata install --user <model name>

To know how to cite a model package::

    macsydata cite <model name>

to have the list of all macsydata subcommand::

    macsydata --help

to have the all options available for a subcommand::

    macsydata <subcommand> --help
