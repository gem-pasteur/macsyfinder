.. _quickstart:


TXSScan Quick Start 
===================

In order to run TXSScan on your favorite dataset as soon as you have installed it, you can simply follow the next steps:

* Type: 
  "``txsscan -h``"
  to see all options available. All command-line options are described in the :ref:`Command-line options section <command-line-label>`.


* On a "metagenomic" dataset for example: 
  "``txsscan all --db_type unordered --sequence_db metagenome.fasta``" 
  to detect all secretion systems and related appendages in a metagenomic dataset.


* On a completely assembled genome (where the gene order is known, and is relevant for systems detection): 
  "``txsscan T2SS Tad --db_type ordered_replicon --sequence_db mygenome.fasta``" 
  to detect the T2SS secretion system and the Tad pilus in a complete genome.

See :ref:`input-dataset-label` for more on input dataset. 


The systems available for detection are the:
    - "Flagellum" -- the bacterial flagellum, involved in motility
    - "T1SS" -- the type 1 secretion system, involved in the secretion of degrading enzymes, toxins,...
    - "T2SS" -- the type 2 secretion system, involved in the secretion of degrading enzymes, toxins,...
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

    Systems have to be spelled as above to run their detection from the command-line. The "all" keyword also allows to detect all available systems in a single run. See the :ref:`Command-line options <command-line-label>`.

    
.. note::

    As some of the above systems have homologous components, the most powerful mode of detection consists in searching "all" the available systems in order to have the most refined assessment, using all the HMM profiles for components detection.


