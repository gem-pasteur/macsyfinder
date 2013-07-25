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


class SystemBank(object):
    """
    build and cached all systems objects. Systems must not be instanciate directly.
    the system_factory must be used. The system factory ensure there is only one instance
    of system for a given name.
    To get a system use the method get_system. if the gene is already cached this instance is returned
    otherwise a new system is build, cached then returned.
    """

    _system_bank = {}


    def __getitem__(self, name):
        """
        :param name: the name of the system
        :type name: string
        :param cfg: the configuration
        :type cfg: :class:`txsscanlib.config.Config` object
        :return: return system corresponding to the name.
         If the system already exists return it otherwise build it an d returni
        :rtype: :class:`txsscanlib.system.System` object
        """
        if name in self._system_bank:
            return self._system_bank[name]
        else:
            raise KeyError(name)


    def __contains__(self, system):
        """
        implement membership test operators
        
        :param system:
        :type system:
        :return: True if the system.name is in , False otherwise
        :rtype: boolean
        """
        return system in self._system_bank.values()

    def __iter__(self):
        """
        """
        return self._system_bank.itervalues()

    def add_system(self, system ):
        """
        :param name: the name of the system
        :type name: string
        :param cfg: the configuration
        :type cfg: :class:`txsscanlib.config.Config` object
        :return: return system corresponding to the name.
         If the system already exists return it otherwise build it an d return
        :rtype: :class:`txsscanlib.system.System` object
        :raise: KeyError if a system with the same name is already registered
        """
        if system in self._system_bank:
            raise KeyError, "a system named %s is already registered" % system.name
        else:
            self._system_bank[system.name] = system

system_bank = SystemBank()


class System(object):
    """
    handle a secretion system.
    """

    def __init__(self, cfg, name, inter_gene_max_space, min_mandatory_genes_required = None, min_genes_required = None):
        """
        :param cfg: the configuration
        :type cfg: :class:`txsscanlib.config.Config` object
        :param name: the name of the system
        :type name: string
        :param inter_gene_max_space: the maximum distance between two genes
        :type inter_gene_max_space: integer
        :param min_mandatory_genes_required: the quorum of mandatory genes to defined this system
        :type min_mandatory_genes_required: integer
        :param min_genes_required: the quorum of genes to defined system
        :type min_genes_required: integer
        """
        self.cfg = cfg
        self.name = name
        self._inter_gene_max_space = inter_gene_max_space
        self._min_mandatory_genes_required = min_mandatory_genes_required
        self._min_genes_required = min_genes_required
        if self._min_mandatory_genes_required is not None and self._min_genes_required is not None:
            if self._min_genes_required < self._min_mandatory_genes_required:
                raise ValueError("min_genes_required must be greater or equal than min_mandatory_genes_required")
        self._mandatory_genes = []
        self._allowed_genes = []
        self._forbidden_genes = []

    @property
    def inter_gene_max_space(self):
        """
        :return: set the maximum distance allowed between 2 genes for this system
        :rtype: integer
        """
        cfg_inter_gene_max_space = self.cfg.inter_gene_max_space(self.name)
        if cfg_inter_gene_max_space is not None:
            return cfg_inter_gene_max_space
        return self._inter_gene_max_space

    @property
    def min_mandatory_genes_required(self):
        """
        :return: get the quorum of mandatory genes required for this system
        :rtype: integer
        """
        cfg_min_mandatory_genes_required = self.cfg.min_mandatory_genes_required(self.name)
        if cfg_min_mandatory_genes_required is not None:
            return cfg_min_mandatory_genes_required
        elif self._min_mandatory_genes_required is None:
            return len(self._mandatory_genes)
        else:
            return self._min_mandatory_genes_required


    @property
    def min_genes_required(self):
        """
        :return: get the quorum of genes required for this system
        :rtype: integer
        """
        cfg_min_genes_required = self.cfg.min_genes_required(self.name)
        if cfg_min_genes_required is not None:
            return cfg_min_genes_required
        elif self._min_genes_required is None:
            return len(self._mandatory_genes)
        else:
            return self._min_genes_required

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

      
