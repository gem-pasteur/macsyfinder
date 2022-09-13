.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2022 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  

.. _installation:


************
Installation
************

MacSyFinder works with models for macromolecular systems that are not shipped with it, 
you have to install them separately. See the :ref:`macsydata section <macsydata>` below.
We also provide container so you can use macsyfinder directly.

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
 
 
MacSyFinder also relies on six Python library dependencies:

 - colorlog
 - colorama
 - pyyaml
 - packaging
 - networkx
 - pandas

These dependencies will be automatically retrieved and installed when using `pip` for installation (see below). 
 

==================================
MacSyFinder Installation procedure
==================================

It is recommended to use `pip` to install the MacSyFinder package.

Archive overview
================

* **doc** => the documentation in html and pdf
* **test** => all what is needed for unitary tests
* **macsypy** => the macsyfinder python library
* **setup.py** => the installation script
* **setup.cfg** => the installation script
* **pyproject.toml** => the project installation build tool
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

.. code-block:: text

    python3 -m pip install macsyfinder


If you do not have the privileges to perform a system-wide installation,
you can either install it in your home directory or
use a `virtual environment <https://virtualenv.pypa.io/en/stable/>`_.

installation in your home directory
"""""""""""""""""""""""""""""""""""

.. code-block:: text

    python3 -m pip install --user macsyfinder


installation in a virtualenv
""""""""""""""""""""""""""""

.. code-block:: text

    python3 -m venv macsyfinder
    cd macsyfinder
    source bin/activate
    python3 -m pip install macsyfinder

To exit the virtualenv just execute the `deactivate` command.
To run `macsyfinder`, you need to activate the virtualenv:

.. code-block:: text

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

To uninstall MacSyFinder (the last version installed), run

.. code-block:: text

  (sudo) pip uninstall macsyfinder

If you install it in a virtualenv, just delete the virtual environment.
For instance if you create a virtualenv name macsyfinder

.. code-block:: text

    python3 -m venv macsyfinder

To delete it, remove the directory

.. code-block:: text

    rm -R macsyfinder

From Conda/Mamba
================

From version 2.0, MacSyFinder is packaged for Conda/Mamba

.. code-block:: text

    mamba install -c macsyfinder=x.x

Where `x.x` is the macsyfinder version you want to install

From container
==============

With Docker
-----------

The docker image is available on Docker Hub (https://hub.docker.com/repository/docker/gempasteur/macsyfinder)
The computations are performed under msf user in /home/msf inside the container.
So You have to mount a directory from the host in the container to exchange data (inputs data, and results) from the host and the container.
The shared directory must be writable by the *msf* user or overwrite the user in the container by your id (see example below)

Furthermore the models are no longer packaged along macsyfinder.
So you have to install them by yourself.
For that we provide a command line tool macsydata which is inspired by pip.

.. code-block:: text

    macsydata search PACKNAME
    macsydata install PACKNAME== or >=, or ... VERSION

To work with Docker you have to install models in a directory which will be mounted in the image at run time

.. code-block:: shell

    mkdir shared_dir
    cd shared_dir

install desired models in my_models directory

.. code-block:: shell

    docker run -v ${PWD}/:/home/msf -u $(id -u ${USER}):$(id -g ${USER})  gempasteur/macsyfinder:<tag> macsydata install --target /home/msf/my_models <MODELS_PACK>

run msf against all models contains in <MODELS_PACK>

.. code-block:: shell

    docker run -v ${PWD}/:/home/msf -u $(id -u ${USER}):$(id -g ${USER})  gempasteur/macsyfinder:<tag> macsyfinder --db-type unordered_replicon --models-dir=/home/msf/my_models/ --models  <MODELS_PACK>  all --sequence-db my_genome.fasta -w 12



With Apptainer (formely Singularity)
------------------------------------

As the docker image is registered in docker hub you can also use it directly with Apptainer (https://apptainer.org/).
Unlike docker you have not to worry about shared directory, your HOME and /tmp are automatically shared.

.. code-block:: shell

    # install desired models in my_models directory
    apptainer run -H ${HOME} docker://gempasteur/macsyfinder:<tag> macsydata install --target my_models <MODELS_PACK>

    # run msf against all models contains in <MODELS_PACK>
    apptainer run -H ${HOME} docker://gempasteur/macsyfinder:<tag> macsyfinder --db-type unordered_replicon --models-dir=my_models --models <MODELS_PACK> all --sequence-db my_genome.fasta -w 12

If you intend to run *apptainer* from host which cannot access internet (cluster node for instance),
you have to

#. download the image locally
#. transfert the image file on the right file system
#. and then use it.

.. code-block:: shell

    apptainer build msf-<tag>.simg docker://gempasteur/macsyfinder:<tag>
    cp msf-<tag>.simg <cluster_file_system>
    apptainer run -H ${HOME} msf-<tag>.simg macsyfinder --db-type unordered_replicon --models-dir=my_models --models <MODELS_PACK> all --sequence-db my_genome.fasta -w 12


.. _macsydata:

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

To list all `macsydata` subcommands::

    macsydata --help

To list all available options for a subcommand::

    macsydata <subcommand> --help

For models not stored in *macsy-models* the commands *available*, *search*,
*installation* from remote or *upgrade* from remote are **NOT** available.

For models **NOT** stored in *macsy-models*, you have to manage them semi-manually.
Download the archive (do not unarchive it), then use *macsydata* to install the archive.