.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _quickstart:


MacSyFinder Quick Start 
=======================

.. only:: html

    .. figure:: ../_static/under_construction.gif

        This page is still under construction

.. only:: latex

    .. figure:: ../_static/under_construction.jpeg

        This page is still under construction


1. We recommend to install MacSyFinder using `pip` in a virtual environment (for further details see :ref:`installation`).

   .. code-block:: bash

        python3 -m venv MacSyFinder
        cd MacSyFinder
        source bin/activate
        pip install macsyfinder

   .. warning::

        `hmmsearch` from the HMMER package (http://hmmer.org/) must be installed.

2. Prepare your data. You need a file containing all protein sequences of your genome of interest in fasta format
   (for further details see :ref:`input-dataset-label`).

3. You need to have models to search in your input data.
   Please refer to :ref:`model_definition` to create your own package of models.
   We will soon provide a set of predefined models for you to test.

4. Command lines:

    * Type:
      :code:`macsyfinder -h`

      To see all options available. All command-line options are described in the :ref:`Command-line options section <command-line-label>`.
      In order to run MacSyFinder on your favorite dataset as soon as you have installed it, you can simply follow the following steps:

    * On a "metagenomic" (unordered) dataset for example:
      :code:`macsyfinder --db-type unordered --sequence-db metagenome.fasta --models-dir my-models --models model_family all`

      will detect all models of model_family modelled in .xml files placed in the *"my-models"* folder without taking into account any gene order.


    * On a completely assembled genome (where the gene order is known, and is relevant for systems' detection):

      :code:`macsyfinder --db-type ordered-replicon --sequence-db mygenome.fasta --models-dir my-models --models model_family ModelA ModelB`

      will detect the macromolecular systems described in the two models *"ModelA"* and *"ModelB"*
      in a complete genome from the *"ModelA.xml"* and *"ModelB.xml"*
      definition files placed in the folder *"my-models/model_family/definitions"*.

.. note::

    Systems names have to be spelled in a case-sensitive way to run their detection from the command-line.
    The name of the System corresponds to the suffix defined for xml files (.xml by default),
    for example *"toto"* for a model defined in *"toto.xml"*.
    
    The *"all"* keyword allows to detect all models available in the definitions folder in a single run.
    See the :ref:`Command-line options <command-line-label>`.



