.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014  Institut Pasteur, Paris.                           
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
    Macsyfinder documentation master file, created by sphinx-quickstart 


Welcome to MacSyFinder's documentation! 
=======================================


	
  MacSyFinder is a program to model and detect macromolecular systems, genetic pathways... in protein datasets. In prokaryotes, these systems have often evolutionarily conserved properties: they are made of conserved components, and are encoded in compact loci (conserved genetic architecture). The user models these systems with MacSyFinder to reflect these conserved features, and to allow their efficient detection. 
  
  Criteria for systems detection include **component content (quorum)**, and **genomic co-localization**. Each component corresponds to a hidden Markov model (HMM) protein profile to perform homology searches with the program Hmmer. 
   
  In order to model macromolecular systems, the user:
    - builds or gather from databanks **HMM protein profiles** for components of interest, 
    - defines **decision rules** for each system in a dedicated XML grammar (see :ref:`system_definition`).


    .. image:: images/figure_main-2.*
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
   macsyview
   
MacSyFinder functioning    
======================= 
.. toctree::
   :maxdepth: 2

   system_definition
   implementation 
   functioning


MacSyFinder API documentation    
=============================
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
   macsypy_error

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

