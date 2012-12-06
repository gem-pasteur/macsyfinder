# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 29, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import os
import logging
_log = logging.getLogger('txsscan.' + __name__)

from subprocess import Popen

class Gene(object):
    """
    handle Gene of a secretion system
    """


    def __init__(self, name, cfg ):
        """
        handle gene
        
        :param name: the name of the gene
        :type name: string
        :param cfg: the configuration 
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.name = name
        self.profile = Profile(self, cfg)# le nom du profile n'est pas deductible du nom de gene?
        self.homologs = []
        self._system = None
        
    def __str__(self):
        s = "name : %s" % self.name
        if self.homologs:
            s += "\n    homologs: "
            for h in self.homologs:
                s += h.name +", "
            s = s[:-2]
        return s
    
    @property
    def system(self):
        """
        :return: the secretion system to which this gene belongs
        :rtype: :class:`txsscanlib.secretion.System` object
        """
        return self._system
    
    @system.setter
    def system(self, system):
        """
        set the system to which this gene belongs
        :param system: the system to which this gene belongs
        :type system: :class:`txsscanlib.secretion.System` object
        """
        self._system = system
        for gene in self.homologs:
            gene.system = system
    
    
    def add_homolog(self, gene):
        """
        add a homolog gene
        
        :param gene: homolg gene to add
        :type gene:  :class:`txsscanlib.secretion.Gene` object
        """
        self.homologs.append(gene)
    

    def get_homologs(self):
        """
        :return: The homolgs genes
        :rtype: list of :class:`txsscanlib.secretion.Gene` object
        """
        return self.homologs

    
class Profile(object):
    """
    handle profile
    """
    
    def __init__(self, gene , cfg):
        """
        
        :param gene: the gene corresponding to this profile
        :type gene_name: :class:`txsscanlib.secretion.Gene` object
        :param cfg: the configuration 
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.gene = gene 
        path = os.path.join(cfg.profile_dir , self.gene.name + cfg.profile_suffix)
        if not os.path.exists(path):
            raise Exception( "%s: No such profile" % path)
        self.path = path
        self.len = self._len()
        self.cfg = cfg 
    
    
    def __len__(self):
        """
        :return: the length of the HMM profile
        :rtype: int
        """
        return self.len
    
    def _len(self):
        """
        Parse the HMM profile to get the length and cache it
        this private method is called at the Profile init
        """
        with open(self.path) as f:
            for l in f:
                if l.startswith("LENG"):
                    length = int(l.split()[1])
                    break
        return length
    
    def __str__(self):
        return "%s : %s" % (self.gene.name, self.path)


    def execute(self):
        """
        execute the hmmsearch with this profile
        :return: ?
        :rtype: ? 
        """
        output_file = os.path.join( self.cfg.res_search_dir, self.gene.name + self.cfg.res_search_suffix )
        err_file = os.path.join( self.cfg.res_search_dir, self.gene.name + os.path.splitext(self.cfg.res_search_suffix)[0]+".err" )
       
        options = { "hmmer_exe" : self.cfg.hmmer_exe,
                   "output_file" : output_file ,
                   "e_value_res" : self.cfg.e_value_res,
                   "profile" : self.gene.name,
                   "sequence_db" : self.cfg.sequence_db,
                   }
        command = "%(hmmer_exe)s -o %(output_file)s -E %(e_value_res)d %(profile)s %(sequence_db)s" % options
        _log.debug( "%s hmmer command line : %s" % (self.gene.name, command) )
        try:
            hmmer = Popen( command ,
                           shell = True ,
                           stdout = None ,
                           stdin  = None ,
                           stderr = err_file ,
                           close_fds = True ,
                           )
        except OSError, err:
            msg = "hmmer execution failed: command = %s : %s" % ( command , err)
            _log.critical( msg, exc_info = True )
            raise err
            
        hmmer.wait()
        if hmmer.returncode != 0:
            msg = "an error occurred during hmmer execution: command = %s : return code = %d check %s" % (command, hmmer.returncode, err_file)
            _log.critical( msg, exc_info = True )
            raise RuntimeError(msg)
        
        return HMMReport(self.gene, output_file, err_file, hmmer.returncode, self.cfg )
    

class HMMReport(object):
    """
    handle HMM report. extract a synthetic report from the raw hmmer output
    """
    
    def __init__(self, gene, hmmer_output, hmmer_err, hmmer_returncode, cfg ):
        """
        :param gene: the gene corresponding to this profile
        :type gene: :class:`txsscanlib.secretion.Gene` object
        :param hmmer_output: The path to hmmer output file
        :type hmmer_output: string
        :param hmmer_err: The path to the hmmer error file
        :type hmmer_err: string
        :param hmmer_returncode: the return code of the hmmer 
        :type hmmer_returncode: int
        :param cfg: the configuration 
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.gene = gene
        self._hmmer_raw_out = hmmer_output
        self._hmmer_err = hmmer_err
        self._hmm_returncode = hmmer_returncode
        self._hmm_extracted = None
        self.cfg = cfg
        
        
    def extract(self):
        """
        Parse the output file of hmmer and produced a new synthetic report file 
        """
        pass
    



class System(object):
    """
    
    """
    
    def __init__(self, name, cfg):
        """
        
        :param name: the name of the system
        :type name: string
        :param cfg: the configuration
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.cfg = cfg
        self.name = name
        self._mandatory_genes = [] #utilser des OrderedDict ?? voir dans le reste du code si on veux y acceder via le nom directement
        self._allowed_genes = []
        self._forbidden_genes = []
    
    def add_mandatory_gene(self, gene):
        """
        add a gene in the list of mandatory genes
        
        :param gene: gene which are mandatory for this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
        gene.system = self
        self._mandatory_genes.append(gene)
    
    def add_allowed_gene(self, gene):
        """
        add a gene in the list of allowed genes
        
        :param gene: gene which should be present in this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
        gene.system = self
        self._allowed_genes.append(gene)
    
    def add_forbidden_gene(self, gene):
        """
        add a gene in the list of forbidden genes
        
        :param gene: gene which must not be present in this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
        gene.system = self
        self._forbidden_genes.append(gene)
    
    @property
    def mandatory_genes(self):
        """
        :return: the list of genes which are mandatory for this secretion system 
        :rtype: list of :class:`txsscanlib.secretion.Gene` object
        """
        return self._mandatory_genes


    @property
    def allowed_genes(self):
        """
        :return: the list of genes which should be present in this secretion system 
        :rtype: list of :class:`txsscanlib.secretion.Gene` object
        """
        return self._allowed_genes
    
    @property
    def forbidden_genes(self):
        """
        :return: the list of genes which cannot be present in this secretion system 
        :rtype: list of :class:`txsscanlib.secretion.Gene` object
        """
        return self._forbidden_genes
    
    