.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2023 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _quickstart:


MacSyFinder Quick Start
=======================
..
    This block is commented (does not appear in compile version)
    .. only:: html

        .. figure:: ../_static/under_construction.gif

            This page is still under construction

    .. only:: latex

        .. figure:: ../_static/under_construction.jpeg

            This page is still under construction


1. We recommend to install MacSyFinder using `pip` in a virtual environment (for further details see :ref:`user_installation`).

   .. code-block:: bash

        python3 -m venv MacSyFinder
        cd MacSyFinder
        source bin/activate
        pip install macsyfinder

   .. warning::

        `hmmsearch` from the HMMER package (http://hmmer.org/) must be installed.

2. Prepare your data. You need a file containing all protein sequences of your genome of interest in a FASTA file
   (for further details see :ref:`input-dataset-label`). In the best case scenario, they would be ordered as the
   corresponding genes are ordered along the replicons.

3. You need to install, or make available to MacSyFinder the models to search in your input genome data.
   Please refer to :ref:`model_definition` to create your own package of models.
   Otherwise, macsy-models contributed by the community are available here: https://github.com/macsy-models
   and can be retrieved and installed using the :ref:`macsydata <macsydata>` command, installed as part of the MacSyFinder suite.


4. Command lines:

    - Type:
      :code:`macsyfinder -h`

      To see all the options available. All command-line options are described in the :ref:`Command-line options section <command-line-label>`.
      In order to run MacSyFinder on your favorite dataset as soon as you have installed the macsy-model of interest,
      you can simply follow the following steps:


    - Install the macsy-models of interest from the `Macsy Models repository <https://github.com/macsy-models>`_:

      :code:`macsydata install some-public-models`


    - On a "unordered" genome dataset:

      :code:`macsyfinder --db-type unordered --sequence-db unordered_genome.fasta --models model_family all`

      will search for systems corresponding to all the models of `model_family` modeled in .xml files shipped with the *"some-public-models"*
      macsy-model package, without taking into account the gene order.

    - On a completely assembled genome (where the gene order is known):

      :code:`macsyfinder --db-type ordered_replicon --sequence-db mygenome.fasta --models-dir my-models --models model_family ModelA ModelB`

      will detect the macromolecular systems described in the two models *"ModelA"* and *"ModelB"*
      in a complete genome from the *"ModelA.xml"* and *"ModelB.xml"*
      definition files placed in the folder *"my-models/model_family/definitions"*.

    - If you want to run the same analysis as above but with local macsy-models not installed by macsydata:

      :code:`macsyfinder --db-type ordered_replicon --sequence-db mygenome.fasta --models-dir my-models --models model_family ModelA ModelB`

      `my-models` is the directory containing  the macsy-model packages.
      NB: The models must follow the :ref:`macsy-models package<package_structure>` structure.

.. note::

    Systems names have to be spelled in a case-sensitive way to run their detection from the command-line.
    The name of the System corresponds to the suffix defined for xml files (.xml by default),
    for example *"toto"* for a model defined in *"toto.xml"*.

    The *"all"* keyword allows to detect all models available in the definitions folder in a single run.
    See the :ref:`Command-line options <command-line-label>`.


An example data set
===================

	We provide `here <https://doi.org/10.6084/m9.figshare.21581280.v1>`_ an example dataset comprising a replicon
	and the output files expected with MacSyFinder, release 2.0 when running the TXSScan macsy-models.
	The genomic dataset consists in the complete sequence of chromosome I from `Vibrio cholerae` O1 biovar El Tor str. N16961
	(published here: https://pubmed.ncbi.nlm.nih.gov/10952301/).

	The chromosome to annotate is presented as a multi-FASTA file of the proteins ordered as the genes encoding them.
	An annotation of the protein secretion systems and appendages was run on the genome, using the macsyfinder set of models ("macsy-model") TXSScan, V1.1.1 in the case of these examples.
	There are two output files offered, the one expected with the "ordered" genome mode of annotation, and the other with the "unordered" mode of genome annotation.
	The following command lines were used to obtain the output files:

	1. The genome is downloaded from `here <https://doi.org/10.6084/m9.figshare.21581280.v1>`_.
	It will serve as an input file in the next command-line examples.

	2. The TXSScan models for annotation of secretion systems are installed.
	The command line is the following:

	:code:`macsydata install TXSScan`
	`# Installs the latest version of TXSScan`

	3. MacSyFinder is run on the genome, here using 8 workers for the HMM search ("-w 8" option):

		- In "ordered" mode:

	:code:`macsyfinder --sequence-db VICH001.B.00001.C001.fasta -o macsyfinder_TXSScan_VICH001_ordered --models TXSScan all --db-type ordered_replicon -w 8`
	`# specified output folder: macsyfinder_TXSScan_VICH001_ordered`


		- In "unordered" mode:

	:code:`macsyfinder --sequence-db VICH001.B.00001.C001.fasta -o macsyfinder_TXSScan_VICH001_unordered --models TXSScan all --db-type unordered -w 8`
	`# specified output folder: macsyfinder_TXSScan_VICH001_unordered`

	The documentation on the generated output files can be consulted :ref:`here <outputs>`.
	See also our FAQ: :ref:`faq-search-mode`

.. note::

	A more comprehensive example of genome datasets with dedicated command lines and expected output files
	can be found `here. <https://doi.org/10.6084/m9.figshare.21716426.v1>`_
