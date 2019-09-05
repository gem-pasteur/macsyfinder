.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014  Institut Pasteur, Paris.                           
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _implementation:

MacSyFinder implementation overview
===================================

MacSyFinder is implemented with an object-oriented architecture. The objects are described in the current :ref:`API documentation <config>`. An overview of the main classes used to model the systems to be detected is provided below.
  
.. digraph:: system_overview

     "System" -> "Gene" -> "Homolog";
     "Gene" -> "System";
     "Homolog" -> "Gene";
     "Gene" -> "Analog";
     "Analog" -> "Gene";
     "Gene" -> "Profile";
     "Gene" -> "HMMReport" -> "Hit";
     "Hit" -> "Gene";
     "Hit" -> "System";
     "Profile" -> "HMMReport"; 
     
The *"System"* class models the systems to detect and contains a list of instances of the *"Gene"* class, which models each component of a given System. The *"Homolog"* and *"Analog"* classes encapsulate a "Gene" and model respectively relationships of homology and analogy between components. 

A *"Gene"* represents a component from the System and refers to an instance of the *"Profile"* object that corresponds to an hidden Markov model protein profile (used for sequence similarity search with the Hmmer program). 

The *"Config"* class (see the :ref:`config`) handles the program parameters, including Hmmer search parameters, and the set of sequences to query (represented by the "Database" object). 

The *"Database"* stores information on the dataset, including necessary information to detect systems in both linear and circular chromosomes (see the :ref:`database`). 

A set of parsers and object factories are used to fill the objects from command-line and input files (*i.e.* the optional configuration file and the XML files describing the systems), and to ensure their uniqueness and integrity. 

Once these objects are initialized and the detection is launched, Hmmer is executed on the sequences of the database (optionally in parallel) with a unique list of profiles corresponding to the systems to detect. Subsequently, Hmmer output files are parsed, and selected hits (given the search parameters provided) are used to fill *"Hits"* objects, which contain information for the detection of the systems. 

During the treatment of the *"Hits"* for *"Systems"* detection, the occurrences of the systems (*"SystemOccurence"* objects) are filled, and the **decision rules** associated with the systems (quorum and co-localization in the case of an "ordered" dataset) are applied. See the following sections for more details on above objects. 



.. _system-implementation:

*****************
The System object
*****************

The :ref:`System object <system>` represents a macromolecular system to detect. 
It is defined *via* a definition file in XML stored in a dedicated location that can be specified *via* the configuration file, or the command-line (`-d` parameter). See :ref:`system-definition-grammar-label` for more details on the XML grammar. 
 
An object :ref:`SystemParser <system_parser>` is used to build a system object from its XML definition file.

A system is named after the file name of its XML definition.
A system has an attribute `inter_gene_max_space` which is an integer,
and three kind of components are listed in function of their presence in the system:

* The genes that must be present in the genome to define this system ("mandatory").
* The genes that can be present, but do not have to be found in every case ("accessory").
* The genes that must not be present in the system ("forbidden").

.. note:: 
    
    a complete description of macromolecular systems modelling is available in the section :ref:`system_definition`


.. _gene-implementation:

***************
The Gene object
***************

The :ref:`Gene object <gene>` represents genes encoding the protein components of a System. 
Each Gene points out its System of origin (:class:`macsypy.system.System`). A Gene must have a correponding HMM protein profile. These profiles are represented by Profile objects (:class:`macsypy.gene.Profile`), and must be named after the gene name. For instance, the gene *gspD* will correspond to the "gspD.hmm" profile file. See :ref:`profile-implementation`). A Gene has several properties described in the :ref:`Gene API <gene>`. 

A Gene may have Homologs or Analogs. An *"Homolog"* (resp. *"Analog"*) object encapsulates a Gene and has a reference to the Gene it is homolog (resp. *"analog"*) to. See the :ref:`Homolog API <homolog-api>` and :ref:`Analog API <analog-api>` for more details. 

.. warning::
    To optimize computation and to avoid concurrency problems when we search several systems, each gene must be instanciated only once, and stored in a *"gene_bank"*.
    gene_bank is a :class:`macsypy.gene.GeneBank` object. 
    The gene_bank and system_bank are filled by the system_parser (:class:`macsypy.system_parser.SystemParser`)


.. _profile-implementation:

******************
The Profile object
******************

Each *"Gene"* component corresponds to a *"Profile"*. The *"Profile"* object is used for the search of the gene with Hmmer. Thus, a *"Profile"* must match a HMM file, which name is based on the profile name. For instance, the *gspG* gene has the corresponding "gspG.hmm" profile file provided at a dedicated location.  


.. _report-implementation:

******************************
Reporting Hmmer search results
******************************

A *"HMMReport"* (:class:`macsypy.report.HMMReport`) object represents the results of a Hmmer program search on the input dataset with a hidden Markov model protein profile. 
This object has methods to extract and build *"Hits"* that are then analyzed for systems assessment. 

It analyses Hmmer raw outputs, and applies filters on the matches (according to :ref:`Hmmer options<hmmer-options>`). See :ref:`hmmer-outputs-label` for details on the resulting output files. 
For profile matches selected with the filtering parameters, *"Hit"* objects are built (see :ref:`the Hit API <hit-label>`). 

