.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  

.. _FAQ:


**************************
Frequently Asked Questions
**************************

How to report an issue?
-----------------------

If you encounter a problem while running MacSyFinder, please submit an issue on the dedicated page of the `GitHub project <https://github.com/gem-pasteur/macsyfinder/issues>`_

To ensure we have all elements to help, please provide: 

- a concise description of the issue
- the expected behavior VS observed one
- the exact command-line used 
- the version of MacSyFinder used
- the exact error message, and if applicable, the `macsyfinder.log` and `macsyfinder.conf` files
- if applicable, an archive (or link to it) with the output files obtained
- if possible, the smallest dataset there is to reproduce the issue
- if applicable, this would also include the macsy-models (XML models plus HMM profiles) used (or precise version of the models if there are publicly available). Same as above, if possible, please provide the smallest set possible of models and HMM profiles. 

All these will definitely help us to help you! ;-) 


What search mode to be used?
----------------------------

Depending on the type of dataset you have, you will have to adapt MacSyFinder's search mode. 

- If you have a fasta file from a complete genome where **proteins are ordered** according to the corresponding genes' order along the replicon, your dataset is entitled to the most powerful search mode: `ordered_replicon` and use the following option `--db-type ordered_replicon`.

- If you have a fasta file of proteins with no sense of the order of the genes along the chromosome(s) or replicon(s), you will have to use the `unordere` search mode with the following option: `--db-type unordered`

- If you have multiple `ordered` replicons to analyse at once, you will have to follow the `Gembase` convention to name the proteins in the fasta file, so that the original replicons can be assessed from their name: :ref:`see here for a description <gembase_convention>` . 
