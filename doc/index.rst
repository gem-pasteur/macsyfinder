.. Txsscan documentation master file, created by
   sphinx-quickstart on Thu Nov 29 18:18:43 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to TXSScan's documentation! 
===================================
	
  TXSScan is a program that detects bacterial protein secretion systems (T1-T6SS) and related appendages (the flagellum, the type  IV pilus, and the Tad pilus). These systems have evolutionarily conserved features: their component content, and their genetic architecture, often in compact loci. We modelled these systems to reflect these features and to allow their efficient detection. 
  
  Criteria for systems detection thus includes **component content (quorum)**, and **genomic co-localization**. Each component corresponds to a hidden Markov model (HMM) protein profile to perform homology searches with the program Hmmer. 
   
  In order to model these protein secretion systems and their related appendages, we:
    - built **HMM protein profiles** for components of interest, 
    - defined **decision rules** for each system in a dedicated XML grammar (see :ref:`system_definition`).


TXSScan documentation contents   
==============================
.. toctree::
   :maxdepth: 2

   installation 
   quickstart
   implementation
   system_definition
   system_parser  
   functioning
   config
   database 
   outputs
   txsscan_error

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

