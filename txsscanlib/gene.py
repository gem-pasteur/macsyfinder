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
from report import OrderedHMMReport, UnOrderedHMMReport

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
        self.profile = Profile(self, cfg)
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
        :rtype: :class:`txsscanlib.system.System` object
        """
        return self._system
    
    @system.setter
    def system(self, system):
        """
        set the system to which this gene belongs
        :param system: the system to which this gene belongs
        :type system: :class:`txsscanlib.system.System` object
        """
        self._system = system
        for gene in self.homologs:
            gene.system = system
    
    
    def add_homolog(self, gene):
        """
        add a homolog gene
        
        :param gene: homolg gene to add
        :type gene:  :class:`txsscanlib.gene.Gene` object
        """
        self.homologs.append(gene)
    

    def get_homologs(self):
        """
        :return: The homolgs genes
        :rtype: list of :class:`txsscanlib.gene.Gene` object
        """
        return self.homologs


class Homolog(Gene):
    """
    handle homologs
    """
    
    def __init__(self, name, cfg, gene_ref, aligned = False ):
        """
        :param name: the name of the gene
        :type name: string
        :param cfg: the configuration 
        :type cfg: :class:`txsscanlib.config.Config` object.
        :param gene_ref: the gene which this one is homolog.
        :type gene_ref: :class:`txsscanlib.gene.Gene` object.
        :param aligned: if True, this gene overlap totally the sequence of the gene reference. Otherwise it overlap partially. 
        :type aligned: boolean
        """
        super(Homolog, self).__init__(name, cfg)
        self.ref = gene_ref
        self.aligned = aligned
    
    def is_aligned(self):
        """
        :return: Tue if this homolog is align to it's gene of reference, False otherwise.
        :rtype: boolean
        """
        return self.aligned
    
    @property
    def gene_ref(self):
        """
        :return: the gene of reference to this homolog
        :rtype: :class:`txsscanlib.gene.Gene` object
        """
        return self.ref
 
 
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
        path = os.path.join(cfg.profile_dir, self.gene.name + cfg.profile_suffix)
        if not os.path.exists(path):
            raise IOError( "%s: No such profile" % path)
        self.path = path
        self.len = self._len()
        self.cfg = cfg 
        self.hmm_raw_output = None
    
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

        :return: an HMM report
        :rtype:  :class:`txsscanlib.report.HMMReport` object
        """
        output_path = os.path.join( self.cfg.res_search_dir, self.gene.name + self.cfg.res_search_suffix )
        err_path = os.path.join( self.cfg.res_search_dir, self.gene.name + os.path.splitext(self.cfg.res_search_suffix)[0]+".err" )

        with  open(err_path, 'w') as err_file:
            options = { "hmmer_exe" : self.cfg.hmmer_exe,
                        "output_file" : output_path ,
                        "e_value_res" : self.cfg.e_value_res,
                        "profile" : self.path,
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
                               close_fds = False ,
                               )
            except Exception, err:
                msg = "hmmer execution failed: command = %s : %s" % ( command , err)
                _log.critical( msg, exc_info = True )
                raise err

            hmmer.wait()
        if hmmer.returncode != 0:
            msg = "an error occurred during hmmer execution: command = %s : return code = %d check %s" % (command, hmmer.returncode, err_path)
            _log.critical( msg, exc_info = True )
            raise RuntimeError(msg)
        self.hmm_raw_output = output_path
        if self.cfg.ordered_db:
            return OrderedHMMReport(self.gene, output_path, self.cfg )
        else:
            return UnOrderedHMMReport(self.gene, output_path, self.cfg )

            