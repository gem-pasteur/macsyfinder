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

import threading

      
class System(object):
    """
    handle a secretion system.
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
        :rtype: list of :class:`txsscanlib.secretion.Gene` objects
        """
        return self._forbidden_genes

    def search_genes(self, genes):
        """
        for each each genes use the profile to perform an HMM and parse the output
        
        :param genes: the genes to search in the genome
        :type genes: list of :class:`txsscanlib.gene.Gene` objects
        """
        # est ce que ca doit rester une methode ou devenir une fonction
        # dans l'etat pas besoin de rester dans system
        # puisqu'on lui pass les genes en arguments
        
        #pour chaque gene   
           ########### bloc a parallelliser  #############
           #recuperer le profil
           #hmmreport = lancer le hmm (execute)
           #extraire le rapport
           ##################### join ####################
        
        def worker(gene):
            profile = gene.profile
            report = profile.execute()
            report.extract()
            report.save_extract()
            
        for g in genes:
            t = threading.Thread(target = worker, args = (g,))
            t.start()
        main_thread = threading.currentThread()    
        for t in threading.enumerate():
            if t is main_thread:
                continue
            t.join()
            