.. _implementation:

TXSScan implementation overview
===============================

TXSScan was implemented using an object-oriented architecture. The objects are described in the current API documentation. 
  
.. digraph:: system_overview

     "System" -> "Gene" -> "Homolog";
     "Gene" -> "System";
     "Homolog" -> "Gene";
     "Gene" -> "Profile";
     "Gene" -> "HMMReport" -> "Hit";
     "Hit" -> "Gene";
     "Hit" -> "System";
     "Profile" -> "HMMReport"; 
     
The *"System"* class models the nanomachines and contains a list of instances of the "Gene" class, which models each component of a given System. The "Homolog" class encapsulates a "Gene" and models relationships of homology between components. 

*"Gene"* represents a component from the System and refers to an instance of the "Profile" object that corresponds to an HMM profile. 

The "Config" class handles the program parameters, including Hmmer search parameters, and the set of sequences to query (represented by the "Database" object). 

The Database stores information on the dataset, including necessary information to detect systems in both linear and circular chromosomes. 

A set of parsers and object factories are used to fill the objects from command-line and input files (*i.e.* the optional configuration file and the XML files describing the systems), and to ensure their uniqueness and integrity. 

Once these objects are initialized and the detection is launched, Hmmer is executed on the sequences of the database (optionally in parallel) with a unique list of profiles corresponding to the systems to detect. Subsequently, Hmmer output files are parsed, and selected hits (given the search parameters provided) are used to fill "Hits" objects, which contain information for the detection of the systems. 

During the treatment of the "Hits" for "Systems" detection, the occurrences of the systems ("SystemOccurence" objects) are filled, and the **decision rules** associated with the systems (quorum and co-localization in the case of an "ordered" dataset) are applied. 
