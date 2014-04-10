.. Txsscan documentation master file, created by
   sphinx-quickstart on Thu Nov 29 18:18:43 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to TXSScan's documentation! 
===================================

	
  TXSScan is a program that detects bacterial protein secretion systems (T1-T6SS) and related appendages (the flagellum, the type  IV pilus, and the Tad pilus). These systems are evolutionarily conserved. They are made of conserved components, and are often encoded  in compact loci (conserved genetic architecture). We modelled these systems to reflect these features and to allow their efficient detection. 
  
  Criteria for systems detection thus includes **component content (quorum)**, and **genomic co-localization**. Each component corresponds to a hidden Markov model (HMM) protein profile to perform homology searches with the program Hmmer. 
   
  In order to model these protein secretion systems and their related appendages, we:
    - built **HMM protein profiles** for components of interest, 
    - defined **decision rules** for each system in a dedicated XML grammar (see :ref:`system_definition`).


    .. image:: images/figure_main-2.png
     :height: 500px
     :align: left


Running MacSyFinder    
===================
.. toctree::
   :maxdepth: 2

   installation 
   quickstart
   input 
   outputs
   
TXSScan functioning    
===================  
.. toctree::
   :maxdepth: 2

   system_definition
   implementation 
   functioning


TXSScan API documentation    
=========================
.. toctree::
   :maxdepth: 2
 
   config
   database
   system
   system_parser
   gene
   profile
   HMMReport
   search_genes
   search_systems
   txsscan_error

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

