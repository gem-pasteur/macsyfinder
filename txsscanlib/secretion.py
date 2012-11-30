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
        self.profile = Profile(name, cfg)# le nom du profile n'est pas deductible du nom de gene?
        self.homologs = []
        
    def __str__(self):
        s = "name : %s" % self.name
        if self.homologs:
            s += "\n    homologs: "
            for h in self.homologs:
                s += h.name +", "
            s = s[:-2]
        return s
    
    
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
    
    def __init__(self, name , cfg):
        """
        
        :param name: the name of the profile
        :type name: string
        :param cfg: the configuration 
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.name = name 
        path = os.path.join(cfg.profile_dir , name + cfg.profile_suffix)
        if not os.path.exists(path):
            raise Exception( "%s: No such profile" % path)
        self.path = path
        self.len = self._len()
        self._cfg = cfg 
    
    
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
        return "%s : %s" % (self.name, self.path)




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
        self._mandatory_genes = [] #utilser des OrderedDict ?? voir dans le reste du code si on vuex y acceder via le nom directement
        self._allowed_genes = []
        self._forbidden_genes = []
    
    def add_mandatory_gene(self, gene):
        """
        add a gene in the list of mandatory genes
        
        :param gene: gene which are mandatory for this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
        self._mandatory_genes.append(gene)
    
    def add_allowed_gene(self, gene):
        """
        add a gene in the list of allowed genes
        
        :param gene: gene which should be present in this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
        self._allowed_genes.append(gene)
    
    def add_forbidden_gene(self, gene):
        """
        add a gene in the list of forbidden genes
        
        :param gene: gene which must not be present in this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
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
    
    