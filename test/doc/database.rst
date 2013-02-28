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
txsscan need also to have the lenght of each sequences and it's position in the database to compute some scores. so it build an index
this file index is store on the same directory as the sequence database. So if this file is present txsscan use it otherwise build this 
index.

The user can force txsscan to rebuild these index with -idx option on the commandline. 
   
     
 
 
 
 
Database API reference
======================
.. automodule:: txsscanlib.database
   :members:
   :private-members:
   :special-members:

