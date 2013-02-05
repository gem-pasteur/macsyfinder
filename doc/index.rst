.. Txsscan documentation master file, created by
   sphinx-quickstart on Thu Nov 29 18:18:43 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Txsscan's documentation!
===================================

Contents:

  bla bla general sur le projet txsscan

Class diagram
=============

The objects Systems, Gene, Profile must be created via theire respective factory. This allow to have only one object of System, Gene or Profile for a given name.
The Homolog objects are composed of a gene and 2 other properties "gene_ref and "aligned". All Gene methods/attributes can be applied to Homolg objects.  
  
.. digraph:: class_diagram

     "System" -> "Gene" -> "Homolog" ;
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
 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

