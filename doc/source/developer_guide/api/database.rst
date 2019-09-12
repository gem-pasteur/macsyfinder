.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.  
    
.. _database:

********
Database
********

The "database" object handles the indexes of the sequence dataset in fasta format,
and other useful information on the input dataset.

MacSyFinder needs to have the length of each sequence and its position in the database
to compute some statistics on Hmmer hits.
Additionally, for ordered datasets ( db_type = 'gembase' or 'ordered_replicon' ),
MacSyFinder builds an internal "database" from these indexes to store information about replicons,
their begin and end positions, and their topology.

The begin and end positions of each replicon are computed from the sequence file,
and the topology from the parsing of the topology file (--topology-file, see :ref:`topology-files`).

Thus it also builds an index (with .idx suffix) that is stored in the same directory as the sequence dataset.
If this file is found in the same folder than the input dataset, MacSyFinder will use it. Otherwise, it will build it.

The user can force MacSyFinder to rebuild these indexes with the "--idx" option on the command-line. 

  
.. _database_api:

database
========
.. automodule:: macsypy.database
   :members:
   :private-members:
   :special-members:


