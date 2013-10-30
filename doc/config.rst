.. _config:

.. _config-definition-label:

*******************
Configuration files
*******************

Options to run TXSScan can be specified in a configuration file. The API described below handles all configuration options for txsscan.
Three locations are parsed to find configuration files: 
 
 * $PREFIX/etc/txsscan/txsscan.conf
 
 * $(HOME)/.txsscan/txsscan.conf
 
 * ./txsscan.conf  
 
 Moreover these three locations options can be passed on the command line.
 
 Each file can defined options, at the end all options are added. if an option is specified several times.
 
 .. note::
    The precedence rules from the less important to the more important are:
 
    $PREFIX/etc/txsscan/txsscan.conf < $(HOME)/.txsscan/txsscan.conf < ./txsscan.conf < "command-line" options
   
    This means that command-line options will always bypass those from the configuration files. In the same flavor, options altering the definition of systems found in the command-line or the configuration file will always overwhelm values from XML definition files.   
 
 The Config object provides some default values and performs some validations of the values, for instance:
 
 
 The configuration files must follow the Python "ini" file syntax:
 
 In TXSScan, five sections are defined:
 
  * **base** : all information related to the protein dataset under study
  
    * *file* : the path to the dataset in Fasta format (*no default value*)
    * *type* : the type of dataset to handle, four types are supported:
       
        * *unordered* : a set of sequences (*e.g.* a metagenomic dataset)
        * *unordered_replicon* : a set of sequences corresponding to a complete replicon (*e.g.* an unassembled complete genome)
        * *ordered_replicon* : a set of sequences corresponding to a complete replicon ordered (*e.g.* an assembled complete genome)
        * *gembase* : a set of multiple ordered replicons.
        
      (*no default value*)
      
    * *replicon_topology* : the topology of the replicon under study. Two topologies are supported: 'linear' and 'circular' (*default* = 'linear')
      This option will be ignored if the dataset type is not ordered (*i.e.* "unordered_replicon" or "unordered").     
  
  * **system**
  
    * *inter_gene_max_space* = list of system name and integer separated by spaces. These values will supersede the values found in the system definition file.  
    * *min_mandatory_genes_required* = list of system name and integer separated by spaces. These values will supersede the values found in the system definition file.
    * *min_genes_required* = list of system name and integer separated by spaces. These values will supersede the values found in the system definition file.
    
  * **hmmer**
    
    * *hmmer_exe* (default= *hmmsearch* )
    * *e_value_res* = (default= *1* )
    * *e_value_sel* = (default= *0.5* )
    * *coverage_profile* = (default= *0.5* )
  
  * **directories**
    
    * *res_search_dir* = (default= *./datatest/res_search* )
    * *res_search_suffix* = (default= *.search_hmm.out* )
    * *profile_dir* = (default= *./profiles* )
    * *profile_suffix* = (default= *.fasta-aln_edit.hmm* )
    * *res_extract_suffix* = (default= *.res_hmm_extract* )
    * *def_dir* = (default= *./DEF/* )
  
  * **general**
  
    * *log_level*: (default= *debug* ) This corresponds to an integer code:
        ========    ========== 
        Level 	    Numeric value
        ========    ==========
        CRITICAL 	50
        ERROR 	    40
        WARNING 	30
        INFO 	    20
        DEBUG 	    10
        NOTSET 	    0
        ========    ==========
    * *log_file* = (default = txsscan.log in directory of the run)
 
  Example of a configuration file::
  
    [base]
    prefix = /path/to/txsscan/home/
    file = %(prefix)s/test/datatest/prru_psae.001.c01.fasta
    type = gembase
    replicon_topology = circular
    
    [system]
    inter_gene_max_space = T2SS 22 Flagellum 44
    min_mandatory_genes_required = T2SS 6 Flagellum 4
    min_genes_required = T2SS 8 Flagellum 4
    
    [hmmer]
    hmmer_exe = hmmsearch
    e_value_res = 1
    i_evalue_sel = 0.5
    coverage_profile = 0.5

    [directories]
    prefix = /path/to/txsscan/home/
    def_dir = %(prefix)s/data/DEF
    res_search_dir = %(prefix)s/test/datatest/res_search/
    res_search_suffix = .search_hmm.out
    profile_dir = %(prefix)s/data/profiles
    profile_suffix = .fasta-aln_edit.hmm
    res_extract_suffix = .res_hmm_extract

   [general]
   log_level = debug
   
     
 
 
Config API reference
====================
.. automodule:: txsscanlib.config
   :members:
   :private-members:
   :special-members:

