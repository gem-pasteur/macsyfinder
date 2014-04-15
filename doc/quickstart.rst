.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014  Institut Pasteur, Paris.                           
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _quickstart:


MacSyFinder Quick Start 
=======================

In order to run MacSyFinder on your favorite dataset as soon as you have installed it, you can simply follow the next steps:

* Type: 
  "``macsyfinder -h``"
  to see all options available. All command-line options are described in the :ref:`Command-line options section <command-line-label>`.


* On a "metagenomic" dataset for example: 

  "``macsyfinder --db_type unordered --sequence_db metagenome.fasta all``" 
  will detect all systems modelled in .xml files placed in the default definition folder in a metagenomic dataset.

  "``macsyfinder --db_type unordered --sequence_db metagenome.fasta -d mydefinitions/ all``" 
  will detect all systems modelled in .xml files placed in the *"mydefinitions"* folder.

* On a completely assembled genome (where the gene order is known, and is relevant for systems detection): 

  "``macsyfinder --db_type ordered_replicon --sequence_db mygenome.fasta -d mydefinitions/ SystemA SystemB``" 
  will detect the systems *"SystemA"* and *"SystemB"* in a complete genome from *"SystemA.xml"* and *"SystemB.xml"* definition files placed in the folder *"mydefinitions"*.

See :ref:`input-dataset-label` for more on input dataset. 


.. The systems available for detection are the:
    - "Flagellum" -- the bacterial flagellum, involved in motility
    - "T1SS" -- the type 1 secretion system, involved in the secretion of degrading enzymes, toxins,...
    - "T2SS" -- the type 2 secretion system, also involved in the secretion of degrading enzymes, toxins,...
    - "T3SS" -- the type 3 secretion, related to the flagellum and dedicated to the secretion into eukaryotic cells
    - "cT4SS" -- the conjugative type 4 secretion system, involved in the transfer of genetic material to other cells
    - "pT4SSi" -- the MPFi-like T4SS dedicated to protein secretion
    - "pT4SSt" -- the MPFt-like T4SS dedicated to protein secretion
    - "T5aSS" -- the "classical" autotransporter 
    - "T5bSS" -- the "two-partner" secretion system
    - "T5cSS" -- the "trimeric" autotransporter
    - "T6SS" -- the type 6 secretion system, involved in protein secretion into bacterial and eukaryotic cells
    - "T4P" -- the type IV pilus, involved in twitching motility, adhesion to cells,...
    - "Tad" -- the Tad pilus, involved in adhesion,...
    

.. note::

    Systems have to be spelled in a case-sensitive way to run their detection from the command-line. The name of the system corresponds to the suffix defined for xml files (.xml by default), for example *"toto"* for a system defined in *"toto.xml"*. 
    
    The *"all"* keyword allows to detect all systems available in the definition folder in a single run. See the :ref:`Command-line options <command-line-label>`.

    


