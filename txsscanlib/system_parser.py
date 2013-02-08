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

import logging
_log = logging.getLogger('txsscan.' + __name__)

from system import system_factory
from gene import gene_factory
from gene import Homolog

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
        path = os.path.join(self.cfg.def_dir, system_name + ".xml")
        if not os.path.exists(path):
            raise Exception("%s: No such sytem definitions" % path)
        system = system_factory.get_system(system_name, self.cfg)
        tree = ET.parse(path)
        root = tree.getroot()
        genes_nodes = root.findall("gene")
        for gene_node in genes_nodes:
            presence = gene_node.get("presence")
            if not presence:
                msg = "Invalid system definition: gene without presence"
                _log.error(msg)
                raise SyntaxError(msg)
            gene = self._parse_gene(gene_node, system)
            if presence == 'mandatory':
                system.add_mandatory_gene(gene)
            elif presence == 'allowed':
                system.add_allowed_gene(gene)
            elif presence == 'forbidden':
                system.add_forbidden_gene(gene)
            else:
                msg = "Invalid system definition: presence value must be either [mandatory, allowed, forbidden] not %s" % presence
                _log.error(msg)
                raise SyntaxError(msg)       
        return system
    
    def _parse_gene(self, node, system):
        """
        parse a xml element gene and build the corresponding object
        
        :param node: a node corresponding to gene element
        :type node: :class:`xml.etree.ElementTree.Element` object 
        :return: the gene object corresponding to the node
        :rtype: :class:`txsscanlib.secretion.System` object 
        """
        
        name = node.get("name")
        if not name:
            msg = "Invalid system definition: gene without name"
            _log.error(msg)
            raise SyntaxError(msg)
        loner = node.get("loner")
        if loner in ("1", "true", "True"):
            loner = True
        elif loner in (None, "0" , "false" , "False"):
            loner = False
            
        exchangeable = node.get("exchangeable")
        if exchangeable in ("1", "true", "True"):
            exchangeable = True
        elif exchangeable in (None, "0" , "false" , "False"):
            exchangeable = False    
        gene = gene_factory.get_gene(self.cfg, name, system, loner, exchangeable)
        
        for homolog_node in node.findall("homologs/gene"):
            gene.add_homolog( self._parse_homolog(homolog_node , gene) )
        return gene
    
    
    def _parse_homolog(self, node, gene_ref):
        """
        parse a xml element gene and build the corresponding object
        
        :param node: a node corresponding to gene element
        :type node: :class:`xml.etree.ElementTree.Element` object 
        :return: the gene object corresponding to the node
        :rtype: :class:`txsscanlib.secretion.System` object 
        """
        
        name = node.get("name")
        if not name:
            msg = "Invalid system definition: gene without name"
            _log.error(msg)
            raise SyntaxError(msg)
        aligned = node.get("aligned")
        if aligned in ("1", "true", "True"):
            aligned = True
        elif aligned in (None, "0" , "false" , "False"):
            aligned = False
        else:
            msg = 'Invalid system definition: invalid value for attribute type for gene %s: %s allowed values are "1","true", "True","0" , "false" , "False" '% (aligned, name)
            _log.error(msg)
            raise SyntaxError(msg)
        system_ref = node.get("system_ref")
        if system_ref is None:
            system = gene_ref.system
        else:
            system = system_factory.get_system(system_ref, self.cfg)
        gene = gene_factory.get_gene(self.cfg, name, system)
        homolog = Homolog(gene, gene_ref, aligned)
        for homolog_node in node.findall("homologs/gene"):
            h2 = self._parse_homolog(homolog_node , gene)
            homolog.add_homolog(homolog)
            homolog.add_homolog( self._parse_homolog(homolog_node , gene) )
        return homolog
    