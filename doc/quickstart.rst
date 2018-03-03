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

  "``macsyfinder --db-type unordered --sequence-db metagenome.fasta all``" 
  will detect all systems modelled in .xml files placed in the default definition folder in a metagenomic dataset.

  "``macsyfinder --db-type unordered --sequence-db metagenome.fasta -d mydefinitions/ all``" 
  will detect all systems modelled in .xml files placed in the *"mydefinitions"* folder.

* On a completely assembled genome (where the gene order is known, and is relevant for systems detection): 

  "``macsyfinder --db-type ordered-replicon --sequence-db mygenome.fasta -d mydefinitions/ SystemA SystemB``" 
  will detect the systems *"SystemA"* and *"SystemB"* in a complete genome from *"SystemA.xml"* and *"SystemB.xml"* definition files placed in the folder *"mydefinitions"*.

See :ref:`input-dataset-label` for more on input datasets. 


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


.. _datatest:   

First trial with a test dataset
*******************************

We included a test dataset in the MacSyFinder package. **By default, it will be installed** in /share/macsyfinder or /usr/share/macsyfinder. But it can be located elsewhere if it was specified during installation.  

This dataset consists in the detection of CRISPR-Cas SubTypes with the definitions in the /share/macsyfinder/DEF folder, using the profiles in the /share/macsyfinder/profiles folder. This classification was previously described in `Makarova et al. 2011 <http://www.ncbi.nlm.nih.gov/pubmed/21552286>`_, and the profiles are from  the `TIGRFAM database <http://www.jcvi.org/cgi-bin/tigrfams/index.cgi>`_ (release 13 of August 15 2012) and some of them were specifically designed for CRISPR-Cas classification (`Haft et. al, 2005 <http://www.ncbi.nlm.nih.gov/pubmed/16292354>`_). The definitions are detailed in the MacSyFinder's paper.

As a sequence dataset, we propose three replicons in /share/macsyfinder/sequence_data/datatest_gembase.fasta: 
    - *Escherichia coli* str. K-12 substr. MG1655 chromosome (ESCO001c01a). Genbank accession number: `NC_000913 <http://www.ncbi.nlm.nih.gov/nuccore/NC_000913>`_.
    - *Haloarcula marismortui* ATCC 43049 plasmid pNG400 (HAMA001p04a). Genbank accession number: `NC_006392 <http://www.ncbi.nlm.nih.gov/nuccore/NC_006392>`_.
    - *Legionella pneumophila* str. Paris, complete genome (LEPN003c01a). Genbank accession number: `NC_006368 <http://www.ncbi.nlm.nih.gov/nuccore/NC_006368>`_.

They were concatenated in a single fasta file, following the "gembase" format proposed :ref:`here <gembase_convention>`, and thus MacSyfinder will treat the three different replicons separately for systems inference. 

To run the detection and classification of all subtypes, type::

    "macsyfinder --db-type gembase --sequence-db 
    /share/macsyfinder/sequence_data/datatest_gembase.fasta all"

To run the detection of the Type-IE subtype only, type::

    "macsyfinder --db-type gembase --sequence-db 
    /share/macsyfinder/sequence_data/datatest_gembase.fasta CAS-TypeIE"

A sample topology file is included /share/macsyfinder/sequence_data/datatest_gembase.topology, and follows the convention in :ref:`here <topology-files>`. It allows to specify a different topology "linear" or "circular" for each replicon in the "gembase" format. Otherwise, by default the topology is set to "circular". It can also be specified in the command-line (see the :ref:`Command-line options <command-line-label>`).

To run the detection using the topology file, type::

    "macsyfinder --db-type gembase --sequence-db 
    /share/macsyfinder/sequence_data/datatest_gembase.fasta 
    --topology-file /share/macsyfinder/sequence_data/datatest_gembase.topology all"

Visualizing expected results with MacSyView
*******************************************

To have an idea of what should be detected with the above test dataset, run :ref:`MacSyView <macsyview>`, the web-browser application for MacSyFinder's results visualization. To do that, open the expected JSON result file with MacSyView: /share/macsyfinder/sequence_data/results.macsyfinder.json. 

A screenshot of MacSyView is included :ref:`here <screenshot>`.


