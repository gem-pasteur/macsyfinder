.. _database:

************
Database API
************

The "database" object handles the indexes of the sequence dataset in fasta format, and other useful information on the input dataset. 

      
TXSScan needs several indexes to run, and speed up the analyses.
 
  * index for hmmsearch (Hmmer program)
  * index for TXSScan

hmmsearch needs to index the sequences to speed up the analyses. The indexes are built by the external tools 
*formatdb* (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/ncbi.tar.gz) or *makeblastdb*. TXSScan tries to find formatdb indexes in the same directory as the sequence file. If the indexes are present 
TXSScan uses these index, otherwise it builds these indexes using formatdb or makeblastdb.

TXSScan needs also to have the length of each sequence and its position in the database to compute some statistics on Hmmer hits. Thus it also builds an index (with .idx suffix) that is stored in the same directory as the sequence dataset. If this file is found in the same folder than the input dataset, TXSScan will use it. Otherwise, it will build it.

The user can force TXSScan to rebuild these indexes with the "--idx" option on the command-line. 
   
Additionally, for ordered datasets ( db_type = 'gembase' or 'ordered_replicon' ), TXSScan builds an internal "database" from these indexes to store information about replicons, their begin and end positions, and their topology.
The begin and end positions of each replicon are computed from the sequence file, and the topology from the parsing of the topology file (--topology-file, see :ref:`topology-files`).
  
  
database API reference
======================
.. automodule:: txsscanlib.database
   :members:
   :private-members:
   :special-members:


