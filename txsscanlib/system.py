# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 29, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import logging
_log = logging.getLogger('txsscan.' + __name__)


class SystemFactory(object):
    """
    build and cached all systems objects. Systems must not be instanciate directly.
    the system_factory must be used. The system factory ensure there is only one instance
    of system for a given name.
    To get a system use the method get_system. if the gene is already cached this instance is returned
    otherwise a new system is build, cached then returned.
    
    """        
    
    system_bank = {}
    
    def get_system(self, name, inter_gene_max_space, cfg):
        """
        :param name: the name of the system
        :type name: string
        :param inter_gene_max_space: the maximum distance allowed between 2 genes
        :type inter_gene_max_space: integer
        :param cfg: the configuration
        :type cfg: :class:`txsscanlib.config.Config` object
        :return: return system corresponding to the name.
        If the system already exists return it otherwise build it an d returni
        :rtype: :class:`txsscanlib.system.System` object
        """
        if name in self.system_bank:
            system =  self.system_bank[name]
        else:
            system = System(name, inter_gene_max_space, cfg)
            self.system_bank[name] = system
        return system
        
system_factory = SystemFactory()
              
class System(object):
    """
    handle a secretion system.
    """

    def __init__(self, name, inter_gene_max_space, cfg):
        """
        :param name: the name of the system
        :type name: string
        :param inter_gene_max_space: the maximum distance allowed between 2 genes
        :type inter_gene_max_space: integer
        :param cfg: the configuration
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.cfg = cfg
        self.name = name
        self._inter_gene_max_space = inter_gene_max_space
        self._mandatory_genes = []
        self._allowed_genes = []
        self._forbidden_genes = []
    
    @property
    def inter_gene_max_space(self):
        return self._inter_gene_max_space

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
        :rtype: list of :class:`txsscanlib.secretion.Gene` objects
        """
        return self._forbidden_genes

      