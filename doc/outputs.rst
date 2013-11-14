.. _outputs:

*************
Output format
*************

TXSScan provides different types of outputs. At each run, TXSScan creates a new folder, which name is based on a fixed prefix and a random suffix, for instance "txsscan-20130128_08-57-46". TXSScan outputs are stored in this run-specific folder. 

.. _hmmer-outputs-label:

1. Hmmer results outputs 
------------------------
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


2. Systems detection results
----------------------------

Different types of tabular outputs are provided. 



3. Logs and configuration files
-------------------------------

A variety of specific output files are built to store info on the TXSScan execution. 

 * txsscan.conf
 
 * txsscan.log
 
 * txsscan.info


