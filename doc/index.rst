.. Txsscan documentation master file, created by
   sphinx-quickstart on Thu Nov 29 18:18:43 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TXSScan's documentation!
===================================

  TXSScan is a program inlucding a framework to model and detect macromolecular systems. It has been applied to the detection of protein secretion systems and related appendages in diderm bacteria (T1SS-T6SS), using Hmmer to perform homologs searches with provided proteic profiles. 
  Criteria for system detection includes *component content (quorum)*, and *genomic co-localization*. These criteria can be tuned by the user.
   
  Its flexible architecture allows the user to describe its own system for detection purpose by:
    - providing its *own profiles* for a new system, 
    - defining its *own decision rules* for system inference in our dedicated XML grammar.

TXSScan implementation overview
===============

The objects Systems, Gene, Profile must be created via their respective factory. This allows to have only one object System, Gene or Profile for a given name.
The Homolog objects are composed of a gene and 2 other properties "gene_ref and "aligned". All Gene methods/attributes can be applied to Homolog objects.  
  
.. digraph:: system_overview

     "System" -> "Gene" -> "Homolog";
     "Gene" -> "Profile";
     "Gene" -> "HMMReport" -> "Hit";
     
txsscanlib API documentation   
============================
.. toctree::
   :maxdepth: 2
   
   system_definition
   config
   system_parser
   system
   gene
   profile
   HMMReport
   search_genes
   search_systems
   database 
   txsscan_error

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

