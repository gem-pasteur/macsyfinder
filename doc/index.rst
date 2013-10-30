.. Txsscan documentation master file, created by
   sphinx-quickstart on Thu Nov 29 18:18:43 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: images/logo_buen.jpg
   :height: 100px
   :width: 200 px
   :align: left


Welcome to TXSScan's documentation! 
===================================
	
  TXSScan is a program that includes a framework to model and detect macromolecular systems. It has been applied to the detection of protein secretion systems and related appendages in diderm bacteria (T1SS-T6SS), using Hmmer to perform homologs searches with provided proteic profiles. 
  Criteria for system detection includes **component content (quorum)**, and **genomic co-localization**. These criteria can be tuned by the user.
   
  Its flexible architecture allows the user to describe its own system for detection purpose by:
    - providing its **own profiles** for a new system, 
    - defining its **own decision rules** for system inference in our dedicated XML grammar (see :ref:`system-definition-grammar-label`).


Table of Contents   
=================
.. toctree::
   :maxdepth: 2

   implementation
   installation 
   system_definition
   config
   database 
   system
   system_parser
   gene
   profile
   HMMReport
   functioning
   search_genes
   search_systems
   txsscan_error

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

