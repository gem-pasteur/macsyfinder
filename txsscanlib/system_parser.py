# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 12, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================

import os
import xml.etree.ElementTree as ET
from secretion import System
from secretion import Gene

class SystemParser(object):
    """
    build a System instance from System definition write in XML and build a 
    """


    def __init__(self, cfg ):
        """
        Constructor

        :param cfg: the configuration of this run
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.cfg = cfg
    
    def parse(self, system_name ):
        """
        parse a system definition in xml format to build the corresponding object
        
        :param system_name: the name of the system to parse
        :type system_name: string
        :return: the system corresponding to the name 
        :rtype: :class:`txsscanlib.secretion.System` object 
        """
        path = os.path.join( self.cfg.def_dir, system_name + ".xml")
        if not os.path.exists(path):
            raise Exception("%s: No such sytem definitions" % path)
        system = System( system_name , self.cfg)
        tree = ET.parse(path)
        root = tree.getroot()
        mandatory_genes_nodes = root.find("mandatory")
        for gene_node in mandatory_genes_nodes.findall("gene"):
            gene = self._parse_gene(gene_node)
            system.add_mandatory_gene(gene)
        allowed_genes_nodes = root.find("allowed")    
        for gene_node in allowed_genes_nodes.findall("gene"):
            gene = self._parse_gene(gene_node)
            system.add_allowed_gene(gene)
        forbidden_gene_nodes = root.find("forbidden")
        for gene_node in forbidden_gene_nodes.findall("gene"):
            gene = self._parse_gene(gene_node)
            system.add_forbidden_gene(gene)
        return system
    
    def _parse_gene(self, node):
        """
        parse a xml element gene and build the corresponding object
        
        :param node: a node corresponding to gene element
        :type node: :class:`xml.etree.ElementTree.Element` object 
        :return: the gene object corresponding to the node
        :rtype: :class:`txsscanlib.secretion.System` object 
        """
        name = node.get("name")
        gene = Gene( name, self.cfg )
        for homolog_node in node.findall("homologs/gene"):
            homolog = self._parse_gene(homolog_node)
            gene.add_homolog(homolog)
        return gene
        