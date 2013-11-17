.. _quickstart:


TXSScan Quick Start 
===================

In order to run TXSScan on your favorite dataset as soon as you have installed it, you can simply follow the next steps:

* Type: 
"``txsscan -h``"
to see full options.

* On a "metagenomic" dataset for example: 
"``txsscan all --db_type unordered --sequence_db metagenome.fasta``" 
to detect all secretion systems and related appendages in a metagenomic dataset

* On a completely assembled genome (where the gene order is known, and is relevant for systems detection): 
"``txsscan T2SS Tad --db_type ordered_replicon --sequence_db mygenome.fasta``" 
to detect the T2SS secretion system and the Tad pilus in a complete genome

See :ref:`database` for more on input dataset. 
