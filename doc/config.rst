.. _config:

******
config
******

Handle the all configuration options for txsscan.
parse 3 location to find configuration files
 
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
    * *ordered* : 
    
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
  
  * **general**
  
    * *log_level* = (default= *debug* )
    * *log_file* = 
 
  example of a file configuration::
  
    [base]
    ordered = yes
    
    [hmmer]
    hmmer_exe = hmmsearch
    e_value_res = 1
    e_value_sel = 0.5
    coverage_profile = 0.5
    
    [directories]
    res_search_dir = ./datatest/res_search
    res_search_suffix = .search_hmm.out
    profile_dir = ./profiles
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

