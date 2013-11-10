.. _database:

*************
Input dataset
*************

The input dataset must be a set of protein sequences in **Fasta format** (see http://en.wikipedia.org/wiki/FASTA_format). The base section in the configuration file (:ref:`config-definition-label`) can be used to specify the path and the type of dataset to deal with.
 
  Four types of protein datasets are supported:
       
        * *unordered* : a set of sequences (*e.g.* a metagenomic dataset)
        * *unordered_replicon* : a set of sequences corresponding to a complete replicon (*e.g.* an unassembled complete genome)
        * *ordered_replicon* : a set of sequences corresponding to a complete replicon ordered (*e.g.* an assembled complete genome)
        * *gembase* : a set of multiple ordered replicons, which format follows the convention we adopted (see :ref:`gembase_convention`).
      
  For "ordered" ("ordered_replicon" or "gembase") datasets only, TXSScan can take into account the shape of the genome: "linear", or "circular". The default is set to "linear". With the "gembase" format, it is possible to specify a topology per replicon with a topology file (see :ref:`gembase_convention`). 




Implementation
==============

The "database" objet handles the indexes of the sequence dataset in fasta format. TXSScan needs several indexes to run, and speed up the analyses.
 
* index for hmmsearch (HMMer program)
* index for TXSScan

hmmsearch needs to index the sequences to speed up the analyses. The indexes are build by the external tools 
*formatdb* (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/ncbi.tar.gz) or *makeblastdb*. txsscan tries to find formatdb index in the same directory as the sequence file. If the indexes are present 
TXSScan uses these index, otherwise it builds these indexes using formatdb or makeblastdb.
TXSScan needs also to have the length of each sequences and its position in the database to compute some statistics on Hmmer hits. Thus it also builds an index (with .idx suffix) that is stored in the same directory as the sequence dataset. If this file is found in the same folder than the input dataset, TXSScan will use it. Otherwise, it will build it.

The user can force TXSScan to rebuild these indexes with the "--idx" option on the command-line. 
   
Additionally, for ordered datasets ( db_type = 'gembase' or 'ordered_replicon' ), TXSScan builds an internal "database" from these indexes to store information about replicons, their begin and end positions, and their topology.
The begin and end positions of each replicon are computed from the sequence file, and the topology from the parsing of the topology file (--topology-file).
this file has the folowing structure:

* a line begining by '#' is ignored 
* one entry per line
* the strucuture is "replicon_identifier" : topology
* topology must be "linear" or "cirular" 
 
example::
 
   #topology file example
   PRRU001c01 : linear
   PSAE001c01 : circular
  
  
database API reference
======================
.. automodule:: txsscanlib.database
   :members:
   :private-members:
   :special-members:

