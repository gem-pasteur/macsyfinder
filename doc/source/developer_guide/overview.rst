.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2019  Institut Pasteur, Paris.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _overview:

MacSyFinder implementation overview
===================================

MacSyFinder is implemented with an object-oriented architecture.
Below a short glossary to fix the vocabulary used in MacSyFinder

.. glossary::
    :sorted:

    Model family
        A set of models, on the same topic.
        it is composed of several definitions which can be sorted in hierachical structure
        and profiles a profile is a hmm profile file.

    Model
        Is a formal description of a macromolecular system.
        Is composed of a definition and a list of profiles.
        at each gene of the Model must correspond a profile

    ModelDefinition

        Is a definition of model, it's serialize as a xml file

    Cluster

        Is a "contiguous" set of hits.

    System

        It's an occurrence of a specific Model on a replicon.
        Basically, it's a cluster or set of clusters which satisfy the Model quorum.

MacSyFinder project structure
-----------------------------

A brief overview of the files and directory constituting the MacSyFinder project

.. glossary::

    doc
        The project is documented using sphinx.
        All sources files needed to generate this documentation is in the directory *doc*

    etc
        This directory contains a template to configure macsyfinder.
        It's allow to set some configuration available for each run and avoid to specify them
        at each run on the command line.
        This file is in *ini* format.

    macsypy
        This the MacSyFinder python library
        Inside macsypy there is a subdirectory *scripts* which are the entry points for
        `macsyfinder` and `macsydata`

    tests
        The code is tested using `unittests`.
        In *tests* the directory *data* contains all data needed to perform the tests.

    utils
        Contains a binary `setsid` needed macsyfinder to parallelize some steps.
        Usually `setsid` is provides by the system, but some macOS version does not provide it.

    CITATION.yml
        A file indicating how to cite macsyfinder in yaml format.

    CONTRIBUTORS
        A file containing the list of code contributors.

    CONTRIBUTING
        A guide on how to contribute to the project.

    COPYRIGHT
        The macsyfinder copyrights.

    COPYING
        The licencing.
        MacSyFinder is released under GPLv3.

    README.md
        Brief information about the project.

    requirements.txt
        The list of python dependencies needed by macsyfinder.
        do not forget to install hmmsearch which is not handle by python packet manager `pip`

    requirements_dev.txt
        The list of extra dependencies needed if you want to contribute to the code.

    setup.py
        The installation recipe.


MacSyFinder architecture overview
---------------------------------

An overview of the main classes.
  
.. figure:: ../_static/macsyfinder_classes.svg

    The macsyfinder classes diagram.
    The classes are not details. only the main attributes allowing us to understand the interaction are mentioned.

    * in green the modules
    * in orange, the concrete class
    * in red the abstract classes
    * in blue the enumeration
    * in purple the dataclass


MacSyFinder functioning overview
--------------------------------
In this section I'll give you an idea of the macsyfinder functioning at very high grain coarse.

As all program the entrypoint is the main function
The goal of `macsyfinder.main` is to parse the command line.
Then to creates a :ref:`config` object and also initialize the logger.
after that it call main_search_systems which contains the macsyfinder logic

The first main_search_systems task is to create models asked by the user on the commandline.
So a DefinitionParser is instantiated.
and the ModelBank and GeneBank are populated

.. note::
    More models than those expressly asked by the user are created.
    macsyfinder parse also models which referred by the asked models trough the analogs, homologs for instance.

Once all models are created, we gather all genes and search them in the replicons.
This step is done in parallel.
The search is done by profil object associated to each gene and rely on the external software hmmsearch.
The parallelization is ensure by search_genes function
The results of this step is a list of hits.

This list is sorted by position and score.
this list is filtered to keep only one hit for each position,
the one with the best score (position is a gene product in a replicon)

For each model asked by the user, we filter the hit list to keep only the hits related to the model.
Those which are link to mandatory, accessory or forbidden gene included the homologs and analogs.

This hits are clustered based on distance constraints describe in the models:

    * **inter_gene_max_space** : the maximum genes allowed between to genes of a system.
    * **lonner** : allow a gene to participate to system even if it does not clusterize with some other genes.

Then we check if each cluster satisfy the quorum described in the model.

    * **min_mandatory_genes** : the minimum of mandatory genes requisite to have a system.
    * **min_genes_required** : the minimum of genes (mandatory + accessory) requisite to have a system.
    * **forbidden_genes** : no forbidden genes may appear in the cluster.

If the model is multi_loci we generate a combination of the clusters and check the quorum for each combination.
If the cluster or combination satisfy the quorum a System is created otherwise a RejectedCluster.

The Systems from the same replicon are sort against their score.




.. _system-implementation:

****************
The Model object
****************

The :ref:`Model object <model>` represents a macromolecular model to detect.
It is defined *via* a definition file in XML stored in a dedicated location that can be specified *via*
the configuration file, or the command-line (`-d` parameter).
See :ref:`model-definition-grammar-label` for more details on the XML grammar.
 
An object :ref:`ModelDEfinitionParser <definition_parser>` is used to build a model object from its XML definition file.

A model is named after the file tree name of its XML definition.
A model has an attribute `inter_gene_max_space` which is an integer,
and three kind of components are listed in function of their presence in the system:

* The genes that must be present in the genome to define this model ("mandatory").
* The genes that can be present, but do not have to be found in every case ("accessory").
* The genes that must not be present in the model ("forbidden").

.. note:: 
    
    a complete description of macromolecular models modelling is available in the section :ref:`model_definition`


.. _gene-implementation:

***************
The Gene object
***************

The :ref:`Gene object <gene>` represents genes encoding the protein components of a Model.
Each Gene points out its Model of origin (:class:`macsypy.model.Model`).
A Gene must have a corresponding HMM protein profile.
These profiles are represented by Profile objects (:class:`macsypy.profile.Profile`),
and must be named after the gene name. For instance, the gene *gspD* will correspond to the "gspD.hmm" profile file.
See :ref:`profile-implementation`). A Gene has several properties described in the :ref:`Gene API <gene>`.

A Gene may have Homologs or Analogs. An *"Homolog"* (resp. *"Analog"*) object encapsulates a Gene and has a reference
to the Gene it is homolog (resp. *"analog"*) to.
See the :ref:`Homolog API <homolog>` and :ref:`Analog API <analog>` for more details.

.. warning::
    To optimize computation and to avoid concurrency problems when we search several systems,
    each gene must be instanciated only once, and stored in a *"gene_bank"*.
    gene_bank is a :class:`macsypy.gene.GeneBank` object. 
    The gene_bank and system_bank are filled by the system_parser (:class:`macsypy.definition_parser.ModelDefinitionParser`)


.. _profile-implementation:

******************
The Profile object
******************

Each *"Gene"* component corresponds to a *"Profile"*.
The *"Profile"* object is used for the search of the gene with Hmmer.
Thus, a *"Profile"* must match a HMM file, which name is based on the profile name.
For instance, the *gspG* gene has the corresponding "gspG.hmm" profile file provided at a dedicated location.


.. _report-implementation:

******************************
Reporting Hmmer search results
******************************

A *"HMMReport"* (:class:`macsypy.report.HMMReport`) object represents the results of a Hmmer program search on
the input dataset with a hidden Markov model protein profile.
This object has methods to extract and build *"Hits"* that are then analyzed for systems assessment. 

It analyses Hmmer raw outputs, and applies filters on the matches (according to :ref:`Hmmer options<hmmer-options>`).
See :ref:`hmmer-outputs-label` for details on the resulting output files.
For profile matches selected with the filtering parameters, *"Hit"* objects are built (see :ref:`the Hit API <hit-label>`). 

.. only:: html

    tests coverage
    --------------

    http://gem.pages.pasteur.fr/MacSyFinder/coverage