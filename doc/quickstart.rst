.. _quickstart:


TXSScan Quick Start 
===================

In order to run TXSScan on your favorite dataset as soon as you installed it, you can simply follow the next steps:

* Type: 
``txsscan -h`` 
to see full options.

* Run: 
``txsscan all --db_type unordered --sequence_db metagenome.fasta`` 
to detect all secretion systems and related appendages in a metagenomic dataset

* Run: 
``txsscan T2SS Tad --db_type ordered_replicon --sequence_db mygenome.fasta`` 
to detect the T2SS secretion system and the Tad pilus in a complete genome
