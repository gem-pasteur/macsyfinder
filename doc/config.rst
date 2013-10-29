.. _config:

*******************
Configuration files
*******************

Handle all the configuration options for txsscan.
Three locations are parsed to find configuration files: 
 
 * $PREFIX/etc/txsscan/txsscan.conf
 
 * $(HOME)/.txsscan/txsscan.conf
 
 * ./txsscan.conf  
 
 Moreover these 3 locations options can be passed on the command line.
 
 Each file can defined options, at the end all options are added. if an option is specified several times.
 The precedence rules from the less important to the more important are:
 $PREFIX/etc/txsscan/txsscan.conf < $(HOME)/.txsscan/txsscan.conf < ./txsscan.conf < command line option
   
 the Config object provide some default values and perform some validations of the values, for instance:
 
 
 the configuration files must follow the python ini file syntax:
 in txsscan 4 sections are defined:
 
  * **base** : all information related to the genome to study
  
    * *file* : the path to the genome (*no default value*)
    * *type* : the type of the base, 4 types are supported:
       
        * *unordered* : a set of sequences (*e.g.* a metagenomic dataset)
        * *unordered_replicon* : a set of sequences corresponding to a complete replicon (*e.g.* an unassembled complete genome)
        * *ordered_replicon* : a set of sequences corresponding to a complete replicon ordered (*e.g.* an assembled complete genome)
        * *gembase* : a set of multiple ordered replicons.
        
      (*no default value*)
      
    * *replicon_topology* : the topology of the replicon, 2 topologies are supported 'linear', 'circular' (*default* = 'linear')
      this option will be ignored if the base type is not ordered (ordered_replicon or unordered).     
  
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
 
  example of a file configuration::
  
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

