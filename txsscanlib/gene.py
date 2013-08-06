# -*- coding: utf-8 -*-

#====================================
# Created on Nov 29, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#====================================


import os
import logging
_log = logging.getLogger('txsscan.' + __name__)

from subprocess import Popen
from threading import Lock
from report import GembaseHMMReport, GeneralHMMReport


class GeneBank(object):
    """
    cached all genes objects. ensure that genes are instanciated only once
    """
    _genes_bank = {}

    def __getitem__(self, name):
        """
        :param name: the name of the gene
        :type name: string
        :param cfg: the configuration
        :type cfg: :class:`txsscanlib.config.Config` object
        :return: return gene corresponding to the name.
          If the gene already exists return it otherwise build it and return it
        :rtype: :class:`txsscanlib.gene.Gene` object
        """
        if name in self._genes_bank:
            return self._genes_bank[name]
        else:
            raise KeyError(name)


    def __contains__(self, gene):
        """
        implement membership test operators
        
        :param gene:
        :type gene:
        :return: True if the gene.name is in , False otherwise
        :rtype: boolean
        
        """
        return gene in self._genes_bank.values()

    def __iter__(self):
        """
        Return an iterator object on the genes contains in the bank
        """
        return self._genes_bank.itervalues()

    def add_gene(self, gene):
        """

        :param name: the name of the gene
        :type name: string
        :param cfg: the configuration
        :type cfg: :class:`txsscanlib.config.Config` object
        :return: return gene corresponding to the name.
          If the gene already exists return it otherwise build it an d return
        :rtype: :class:`txsscanlib.gene.Gene` object
        :raise: KeyError if a gene with the same name is already registered
        """
        if gene in self._genes_bank:
            raise KeyError("a gene named %s is already registered" % gene.name)
        else:
            self._genes_bank[gene.name] = gene

gene_bank = GeneBank()


class Gene(object):
    """
    handle Gene of a secretion system

    """

    def __init__(self, cfg, name,
                 system,
                 loner = False,
                 exchangeable = False,
                 multi_system = False,
                 inter_gene_max_space = None ):
        """
        handle gene

        :param cfg: the configuration. 
        :type cfg: :class:`txsscanlib.config.Config` object
        :param name: the name of the gene.
        :type name: string.
        :param system: the system which belongs this gene/
        :type system: :class:`txsscanlib.system.System` object.
        :param loner: True if a gene can be isolated on the genome, False otherwise.
        :type loner: boolean.
        :param exchangeable: True if this gene can be replaced with one of its homologs whithout any effects on the system, False otherwise.
        :type exchangeable: boolean.
        :param multi_system: True if a gene can belong to different system. 
        :type multi_system: boolean.
        :param inter_gene_max_space: the maximum space between this gene and an other gene of the system.
        :type inter_gene_max_space: integer
        """
        self.name = name 
        self.profile = profile_factory.get_profile(self, cfg)
        """:ivar profile: The profile HMM Profile corresponding to this gene :class:`txsscanlib.gene.Profile` object"""

        self.homologs = []
        self._system = system
        self._loner = loner
        self._exchangeable = exchangeable
        self._multi_system = multi_system
        self._inter_gene_max_space = inter_gene_max_space

    def __str__(self):
        s = "name : %s" % self.name
        if self.homologs:
            s += "\n    homologs: "
            for h in self.homologs:
                s += h.name + ", "
            s = s[:-2]
        return s

    @property
    def system(self):
        """
        :return: the secretion system to which this gene belongs
        :rtype: :class:`txsscanlib.system.System` object
        """
        return self._system

    @property
    def loner(self):
        """
        :return: True if the gene can be isolated on the genome false otherwise
        :rtype: boolean
        """
        return self._loner

    @property
    def exchangeable(self):
        """
        :return: True if this gene can be replaced with one of its homologs whithout any effects on the system, False otherwise.
        :rtype: boolean.
        """
        return self._exchangeable

    @property
    def multi_system(self):
        """
        :return: True if this gene can belong to different systems, False otherwise.
        :rtype: boolean.
        """
        return self._multi_system
    
    @property
    def inter_gene_max_space(self):
        """
        :return: The maximum distance between this gene and an other from the system. 
                 If the value is not set at the gene level return the value set at the system level.
        :rtype: integer.
        """
        if self._inter_gene_max_space is not None:
            return self._inter_gene_max_space
        else:
            return self._system.inter_gene_max_space
    
    def add_homolog(self, homolog):
        """
        add a homolog gene

        :param homolog: homolog to add
        :type homolog:  :class:`txsscanlib.gene.Homolog` object 
        """
        self.homologs.append(homolog)


    def get_homologs(self):
        """
        :return: The homologs genes
        :type: list of :class:`txsscanlib.gene.Gene` object
        """
        return self.homologs
    
    def __eq__(self, gene):
        """
        :return: True if the profiles (genes) are the same, False otherwise.
        :type gene: :class:`txsscanlib.gene.Gene` object.
        :rtype: boolean.
        """
        return self.name == gene.name
        
        
    def is_homolog(self, gene):
        """
        :return: True if the genes are homologs, False otherwise.
        :type gene: :class:`txsscanlib.gene.Gene` object.
        :rtype: boolean.
        """

        if self == gene:
            return True
        else:
            for h in self.homologs:
                if gene == h.gene:
                    return True
               
        return False       
    
    def is_mandatory(self, system):
        if self in system.mandatory_genes:
            return True
        else:
            return False    
    
    def is_allowed(self, system):
        if self in system.allowed_genes:
            return True
        else:
            return False
            
    def is_forbidden(self, system):
        if self in system.forbidden_genes:
            return True
        else:
            return False
    
    def is_authorized(self, system):
        for m in (system.mandatory_genes+system.allowed_genes):
            if self == m:
                return True
            if m.exchangeable and m.is_homolog(self):
                return True
            
        return False
        
class Homolog(object):
    """
    handle homologs
    """

    def __init__(self, gene, gene_ref, aligned = False ):
        """
        :param gene: the gene
        :type gene: :class:`txsscanlib.gene.Gene` object.
        :param gene_ref: the gene to which this one is homolog.
        :type gene_ref: :class:`txsscanlib.gene.Gene` object.
        :param aligned: if True, the profile of this gene overlaps totally the sequence of the reference gene profile. Otherwise, only partial overlapping between the profiles. 
        :type aligned: boolean
        """
        self.gene = gene 
        """:ivar gene: gene """

        self.ref = gene_ref 
        self.aligned = aligned

    def __getattr__(self, name):
        return getattr(self.gene, name)

    def is_aligned(self):
        """
        :return: True if this homolog is aligned to its gene of reference, False otherwise.
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


class ProfileFactory():
    """
    build and cached all Profile objects. Profiles must not be instanciate directly.
    the profile_factory must be used. The profile_factory ensure there is only one instance
    of profile for a given name.
    To get a profile use the method get_profile. If the profile is already cached this instance is returned
    otherwise a new profile is build, cached then returned.

    """
    _profiles = {}

    def get_profile(self, gene, cfg):
        """
        :return: return profile corresponding to the name.
                 If the profile already exists return it otherwise build it and return
        :rtype: :class:`txsscanlib.gene.Profile` object
        """
        if gene.name in self._profiles:
            profile = self._profiles[gene.name]
        else:
            profile = Profile(gene, cfg)
            self._profiles[gene.name] = profile
        return profile 

profile_factory = ProfileFactory()


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
        self._report = None
        self._lock = Lock()


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
        with self._lock:
            # the results of HMM is cached 
            # so HMMsearch is executed only once per run
            # if this method is called several times the first call induce the execution of HMMsearch and generate a report
            # the other calls return directly this report
            if self._report is not None:
                return self._report
            output_path = os.path.join( self.cfg.working_dir, self.gene.name + self.cfg.res_search_suffix )
            err_path = os.path.join( self.cfg.working_dir, self.gene.name + os.path.splitext(self.cfg.res_search_suffix)[0]+".err" )

            with  open(err_path, 'w') as err_file:
                options = { "hmmer_exe" : self.cfg.hmmer_exe,
                            "output_file" : output_path ,
                            "e_value_res" : self.cfg.e_value_res,
                            "profile" : self.path,
                            "sequence_db" : self.cfg.sequence_db,
                           }
                #command = "%(hmmer_exe)s -o %(output_file)s -E %(e_value_res)d %(profile)s %(sequence_db)s" % options
                command = "%(hmmer_exe)s --cpu 0 -o %(output_file)s -E %(e_value_res)d %(profile)s %(sequence_db)s " % options
                _log.info( "%s hmmer command line : %s" % (self.gene.name, command) )
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
            if self.cfg.db_type == 'gembase':
                report = GembaseHMMReport(self.gene, output_path, self.cfg )
            else:
                report = GeneralHMMReport(self.gene, output_path, self.cfg )
            self._report = report
            return report
            
