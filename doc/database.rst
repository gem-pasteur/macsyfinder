.. _database:

********
database
********

Handle the indexes of the sequence base in fasta format. txsscan need several indexes to work and speed up analyses.
 
* index for hmmsearch
* index for txsscan

hmmsearch need to index the sequence to speed up the annalyses. The indexes are build by the external tools 
formatdb (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/ncbi.tar.gz). txsscan try to find formatdb index in the same directory as the sequence database file. If the indexes are present 
txsacn use these index, otherwise build these index using formatdb.
txsscan need also to have the length of each sequences and it's position in the database to compute some scores. So it build an index,
this file index (with .idx siffix) is stored on the same directory as the sequence database. So if this file is present txsscan use it otherwise build this 
index.

The user can force txsscan to rebuild these indexes with -idx option on the commandline. 
   
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

