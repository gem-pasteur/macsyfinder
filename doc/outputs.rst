.. _outputs:

*************
Output format
*************

MacSyFinder provides different types of outputs. At each run, MacSyFinder creates a new folder, whose name is based on a fixed prefix and a random suffix, for instance "txsscan-20130128_08-57-46". MacSyFinder outputs are stored in this run-specific folder. 

.. _hmmer-outputs-label:

Hmmer results outputs 
---------------------
Raw Hmmer outputs are provided, as long with processed tabular outputs that includes hits filtered as specified by the user. For instance, the Hmmer search for SctC homologs with the corresponding profile will produce as a result two files: "sctC.search_hmm.out" and "sctC.res_hmm_extract". 

The processed output "sctC.res_hmm_extract" recalls on the first lines the parameters used for hits filtering and relevant information on the matches, as 
for instance::

    # gene: sctC extract from /Users/bob/txsscan_results/txsscan-20130128_08-57-46/sctC.search_hmm.out hmm output
    # profile length= 544
    # i_evalue threshold= 0.001000
    # coverage threshold= 0.500000
    # hit_id replicon_name position_hit hit_sequence_length gene_name gene_system i_eval score profile_coverage sequence_coverage begin end
    PSAE001c01_006940       PSAE001c01      3450    803     sctC    T3SS    1.1e-41 141.6   0.588235        0.419676        395     731
    PSAE001c01_018920       PSAE001c01      4634    776     sctC    T3SS    9.2e-48 161.7   0.976103        0.724227        35      596
    PSAE001c01_031420       PSAE001c01      5870    658     sctC    T3SS    2.7e-52 176.7   0.963235        0.844985        49      604
    PSAE001c01_051090       PSAE001c01      7801    714     sctC    T3SS    1.9e-46 157.4   0.571691        0.463585        374     704


.. note::
    Each tabular output file contains a header line describing each column in the output.


Systems detection results
-------------------------

Different types of tabular outputs are provided. Headers are provided with the content of the lines in the file.

  * txsscan.tab - for **"ordered" datasets only** (db_type is "ordered_replicon" or "gembase"). It provides a summary of the number of each type of systems that have been detected. 
  
  * txsscan.report - contains all the sequence identifiers of proteins detected as being part of a system, along with statistics on the Hmmer hit, and the status of the component in the system. 
  
  * txsscan.summary - contains a line of information for each detected system.


txsscan.tab
***********
For each replicon, a line gives the occurrences of systems that were asked for detection. For example, if the detection was run for the Flagellum and the T6SS, the output will look like::

  #Replicon Flagellum_single_locus Flagellum_multi_loci T6SS_single_locus T6SS_multi_loci	
  escherichia06 1   1   1   0

which means that this "escherichia06" genome harbors 1 flagellum in a single locus, 1 flagellum scattered on multiple loci, and 1 T6SS in a single locus. 

txsscan.report
**************
Each line corresponds to a "hit" that has been assigned to a detected system. It includes:
    * Hit_Id - the sequence identifier of the hit
    * Replicon_name	- the name of the replicon it belongs to
    * Position - the position of the sequence in the replicon
    * Sequence_length - the length of the sequence
    * Gene - the name of the components matched with the profile
    * Reference_system - the system that includes the component matched
    * Predicted_system - the system assigned
    * System_Id - the unique identifier attributed to the detected system
    * System_status	- the status of the detected system
    * Gene_status - the status of the component in the assigned system's definition 
    * i-evalue - Hmmer statistics, the indepent-evalue
    * Score	- Hmmer score
    * Profile_coverage - the percentage of the profile covered by the alignment with the sequence
    * Sequence_coverage - the percentage of the sequence covered by the alignment with the profile
    * Begin_match - the position in the sequence where the profile match begins
    * End_match - the position in the sequence where the profile match ends

txsscan.summary
***************
Each line corresponds to a system that has been detected. It includes:
    * Replicon_name	- the name of the replicon 
    * System_Id	- the unique identifier attributed to the detected system
    * Reference_system - the type of system detected	
    * System_status	- the status of the system
    * Nb_loci - the number of loci that constitutes the system
    * Nb_Ref_mandatory - the number of mandatory genes in the system definition
    * Nb_Ref_accessory - the number of accessory genes in the system definition
    * Nb_Ref_Genes_detected_NR - the number of different components (accessory+mandatory) in the system 
    * Nb_Genes_with_match - the number of components detected with the profiles in the system
    * System_length	- the full number of components (with match or not) in the locus (or loci) that constitutes the system 
    * Nb_Mandatory_NR - the number of different mandatory components matched  
    * Nb_Accessory_NR - the number of different accessory components matched 
    * Nb_missing_mandatory - the number of mandatory components from the system definition with no match in this system occurrence
    * Nb_missing_accessory - the number of accessory components from the system definition with no match in this system occurrence	
    * List_missing_mandatory - the list of the missing mandatory components
    * List_missing_accessory - the list of the missing accessory components
    * Loci_positions - the sequence position (rank of the fasta sequence in the input sequence file) of the different loci encoding the system 
    * Occur_Mandatory - counts of the mandatory components
    * Occur_Accessory - counts of the accessory components
    * Occur_Forbidden - counts of the forbidden components

txsscan.json
************
This file is used by MacSyView, for graphical output purpose. It must be loaded through MacSyView to graphically visualize detected systems. For more details, see :ref:`MacSyView's description <macsyview>`.


Logs and configuration files
----------------------------

Two specific output files are built to store information on the MacSyFinder execution: 

 * txsscan.conf - contains the configuration information of the run. It is useful to recover the parameters used for the run. 
 
 * txsscan.log - the log file, contains raw information on the run. Please send it to us with any bug report. 
  


