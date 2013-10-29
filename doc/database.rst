.. _database:

*************
Input dataset
*************

The input dataset must be a set of protein sequences in **Fasta format**. The base section in the configuration file (:ref:`config-definition-label`) can be used to specify the path and the type of dataset to deal with.
 
  Four types of protein datasets are supported:
       
        * *unordered* : a set of sequences (*e.g.* a metagenomic dataset)
        * *unordered_replicon* : a set of sequences corresponding to a complete replicon (*e.g.* an unassembled complete genome)
        * *ordered_replicon* : a set of sequences corresponding to a complete replicon ordered (*e.g.* an assembled complete genome)
        * *gembase* : a set of multiple ordered replicons.
      
  For "ordered" ("ordered_replicon" or "gembase) datasets only, TXSScan can take into account the shape of the genome: "linear", or "circular". The default is set to "linear".    




Implementation
==============

The "database" objet handles the indexes of the sequence base in fasta format. txsscan need several indexes to work and speed up analyses.
 
* index for hmmsearch
* index for txsscan

hmmsearch need to index the sequence to speed up the analyses. The indexes are build by the external tools 
formatdb (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/ncbi.tar.gz). txsscan try to find formatdb index in the same directory as the sequence database file. If the indexes are present 
TXSScan uses these index, otherwise build these index using formatdb.
TXSScan needs also to have the length of each sequences and its position in the database to compute some scores. So it builds an index,
this index file (with .idx suffix) is stored on the same directory as the sequence database. So if this file is present txsscan use it otherwise build this 
index.

The user can force txsscan to rebuild these indexes with -idx option on the command-line. 
   
In addition, for ordered database ( db-type = 'gembase' or 'ordered_replicon' ), of these two kind of indexes, txsscan build an internal "database" to store information 
about replicons, the begin, end and the topology.
The begin and end of each replicon are computed from the sequences base, the topology come from the parsing of the topology file (--topology-file).
this file has the folowing structure:

* a line begining by '#' is ignore 
* one entry per line
* replicon identifier : topology
* topology must be linear or cirular 
 
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

