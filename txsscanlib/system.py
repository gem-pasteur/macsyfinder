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
from txsscan_error import SystemInconsistencyError

class SystemBank(object):
    """
    Build and store all Systems objects. Systems must not be instanciated directly.
    This system factory must be used. It ensures there is a unique instance
    of a system for a given system name.
    To get a system, use the method __getitem__ via the "[]". If the System is already cached in the SystemBank, it is returned.
    Otherwise a new system is built, stored and then returned.
    """

    _system_bank = {}


    def __getitem__(self, name):
        """
        :param name: the name of the system
        :type name: string
        :param cfg: the configuration object
        :type cfg: :class:`txsscanlib.config.Config` object
        :return: the system corresponding to the name.
         If the system already exists, return it, otherwise build it and return it.
        :rtype: :class:`txsscanlib.system.System` object
        """
        if name in self._system_bank:
            return self._system_bank[name]
        else:
            raise KeyError(name)


    def __contains__(self, system):
        """
        Implement the membership test operator
        
        :param system: the system to test
        :type system: :class:`txsscanlib.system.System` object
        :return: True if the system is in the System factory, False otherwise
        :rtype: boolean
        """
        return system in self._system_bank.values()

    def __iter__(self):
        """
        Return an iterator object on the systems contained in the bank
        """
        return self._system_bank.itervalues()
    
    def __len__(self):
        """
        :return: the number of systems stored in the bank
        :rtype: integer
        """
        return len(self._system_bank)
    
    
    def add_system(self, system ):
        """
        :param name: the name of the system
        :type name: string
        :param cfg: the configuration object
        :type cfg: :class:`txsscanlib.config.Config` object
        :return: the system corresponding to the system's name passed as an argument
         If the system already exists, return it. Otherwise, build it and return it. 
        :rtype: :class:`txsscanlib.system.System` object
        :raise: KeyError if a system with the same name is already registered.
        """
        if system in self._system_bank:
            raise KeyError, "a system named %s is already registered in the systems' bank" % system.name
        else:
            self._system_bank[system.name] = system

system_bank = SystemBank()


class System(object):
    """
    Handle a secretion system.
    """

    def __init__(self, cfg, name, inter_gene_max_space, min_mandatory_genes_required = None, min_genes_required = None, max_nb_genes = None, multi_loci = False):
        """
        :param cfg: the configuration object
        :type cfg: :class:`txsscanlib.config.Config` object
        :param name: the name of the system
        :type name: string
        :param inter_gene_max_space: the maximum distance between two genes (**co-localization** parameter)
        :type inter_gene_max_space: integer
        :param min_mandatory_genes_required: the quorum of mandatory genes to define this system
        :type min_mandatory_genes_required: integer
        :param min_genes_required: the quorum of genes to define this system
        :type min_genes_required: integer
        :param max_nb_genes: 
        :type max_nb_genes: integer
        :param multi_loci: 
        :type multi_loci: boolean
        """
        self.cfg = cfg
        self.name = name
        self._inter_gene_max_space = inter_gene_max_space
        self._min_mandatory_genes_required = min_mandatory_genes_required
        self._min_genes_required = min_genes_required
        if self._min_mandatory_genes_required is not None and self._min_genes_required is not None:
            if self._min_genes_required < self._min_mandatory_genes_required:
                raise SystemInconsistencyError("min_genes_required must be greater or equal than min_mandatory_genes_required")
        self._max_nb_genes = max_nb_genes
        self._multi_loci = multi_loci
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
        :return: get the minimum number of genes to assess for the system presence.
        :rtype: integer
        """
        cfg_min_genes_requireds = self.cfg.min_genes_required(self.name)
        if cfg_min_genes_requireds is not None:
            return cfg_min_genes_requireds
        elif self._min_genes_required is None:
            return len(self._mandatory_genes)
        else:
            return self._min_genes_required

    @property
    def max_nb_genes(self):
        """
        :return: the maximum number of genes to assess the system presence.
        :rtype: int (or None)
        """
        cfg_max_nb_genes = self.cfg.max_nb_genes(self.name)
        if cfg_max_nb_genes:
            return cfg_max_nb_genes
        else:
            return self._max_nb_genes

    @property
    def multi_loci(self):
        """
        :return: True if the system is authorized to be inferred from multiple loci, False otherwise
        :rtype: boolean
        """
        cfg_multi_loci = self.cfg.multi_loci(self.name)
        if cfg_multi_loci:
            return cfg_multi_loci
        else:
            return self._multi_loci
        
    def add_mandatory_gene(self, gene):
        """
        Add a gene to the list of mandatory genes

        :param gene: gene that is mandatory for this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
        self._mandatory_genes.append(gene)

    def add_allowed_gene(self, gene):
        """
        Add a gene to the list of allowed genes

        :param gene: gene that is allowed to be present in this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
        self._allowed_genes.append(gene)

    def add_forbidden_gene(self, gene):
        """
        Add a gene to the list of forbidden genes

        :param gene: gene that must not be found in this system
        :type gene: :class:`txsscanlib.secretion.Gene` object
        """
        self._forbidden_genes.append(gene)

    @property
    def mandatory_genes(self):
        """
        :return: the list of genes that are mandatory in this secretion system 
        :rtype: list of :class:`txsscanlib.secretion.Gene` objects
        """
        return self._mandatory_genes

    @property
    def allowed_genes(self):
        """
        :return: the list of genes that are allowed in this secretion system 
        :rtype: list of :class:`txsscanlib.secretion.Gene` objects
        """
        return self._allowed_genes

    @property
    def forbidden_genes(self):
        """
        :return: the list of genes that are forbidden in this secretion system 
        :rtype: list of :class:`txsscanlib.secretion.Gene` objects
        """
        return self._forbidden_genes

      
