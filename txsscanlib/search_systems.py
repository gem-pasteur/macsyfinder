# -*- coding: utf-8 -*-

###############################
# Created on Feb 14, 2013
#
# @author: sabby
# @contact: sabby@pasteur.fr
# @organization: Institut Pasteur
# @license: GPLv3
################################

import threading
import logging
from report import Hit # required? 
#from operator import attrgetter # To be used by the groupby()
#from itertools import groupby
import itertools, operator
_log = logging.getLogger('txsscan.' + __name__)


"""
def search_systems(systems, hits, cfg):
   
    #gets all hits obtained from HMMER runs and applies decision rules to assess system occurences. 
    #Criteria are the quorum of genes, and the colocalization in the case of \"ordered\" datasets. 
    
    
    #_log.debug("Starting system detection with search_systems")
    
    # Deals with different dataset types using Pipeline ?? 
       
    clusters=ClustersHandler()
    
    # Browse positions instead of Hit !!! 
    for hit in sort(hits):
        if hit.is_mandatory() or hit.is_allowed(): # in any system
            cur_system = hit.get_system()
            # will take the appropriate inter_gene_max_space and all hits within this area
            # replies True only if neighbors are allowed in the current system !
            neighbies=hit.get_neighbors(cur_system, cfg)
            #if hit.has_neighbors(cur_system, cfg): 
            if neighbies:
                #hit.get_neighbors(cur_system)
                cluster = Cluster(hit, cfg)
                for n in neighbies:
                    cluster.add(n)
                clusters.add(cluster)
                # Computes where to restart to check
                jump_pos=cluster.end - pos + 1
            else:
                if hit.is_loner():
                    clusters.add_hit(hit)
                    detected_syst.add_gene(cur_system, hit)
                else:
                    pass 
"""




class ClustersHandler(object):
    """
    Deals with sets of clusters found in a dataset. Conceived to store only clusters for a same replicon.
    """
    
    def __init__(self, cfg):
        """
        :param cfg: The configuration object built from default and user parameters.
        :type cfg: :class:`txsscanlib.config.Config` 
        """
        self.clusters = []
        self.replicon_name = ""
    
    def add(self, cluster):
        if not self.replicon_name:
            self.replicon_name = cluster.replicon_name            
            cluster.save()
            self.clusters.append(cluster)
        elif self.replicon_name == cluster.replicon_name:            
            cluster.save()
            self.clusters.append(cluster)
        else:
            msg = "Attempting to add a cluster in a ClustersHandler dedicated to another replicon !"
            _log.critical(msg)
            raise Exception(msg)
        
    def __str__(self):
        to_print=""
        for cluster in self.clusters:
            to_print+=str(cluster)
            
        return to_print

                        
class Cluster(object):
    """
    Stores a set of contiguous hits. 
    """
    
    def __init__(self, cfg):
        self.hits = []
        self.systems = []
        self.replicon_name = ""
        self.begin = 0
        self.end = 0
        self.state = ""
        self.putative_system = ""
        #self.add(hit)
    
    def __len__(self):
        return len(self.hits)
    
    def __str__(self):
        
        pos = []
        seq_ids = []
        gene_names = [] 
        for h in self.hits:
            pos.append(h.position)
            seq_ids.append(h.id)
            gene_names.append(h.gene.name)
            
        #return "--- Cluster ---\n%s\n%s\n%s"%(str(seq_ids), str(gene_names), str(pos))
        if self.state == "clear":
            return "--- Cluster %s ---\n%s\n%s\n%s"%(self.putative_system, str(seq_ids), str(gene_names), str(pos))
        else:
            return "--- Cluster %s ? ---\n%s\n%s\n%s"%(self.putative_system, str(seq_ids), str(gene_names), str(pos))
        
    def add(self, hit):
        # need to update cluster bounds
        if len(self.hits) == 0:
            self.begin = hit.get_position()
            self.end = self.begin
            self.replicon_name = hit.replicon_name
            self.hits.append(hit)
        else:
            if(self.replicon_name == hit.replicon_name):
                if hit.get_position() < self.begin:
                    self.begin = hit.get_position()
                elif hit.get_position() > self.end:
                    self.end = hit.get_position()
                else:
                    if hit.get_position() > self.begin and hit.get_position() < self.end:
                        _log.debug("Weird cluster inclusion hit : %s"%hit)
                            
                self.hits.append(hit)
                
            else:
                msg = "Attempting to gather in a cluster hits from different replicons ! "
                _log.critical(msg)
                raise Exception(msg)
                    
    def save(self):
        
        if not self.putative_system :
            systems={}
            for h in self.hits:
                syst=h.system.name
                if not systems.has_key(syst):
                    systems[syst]=1
                else:
                    systems[syst]+=1
            
            #print systems
                
            max_syst=0
            tmp_syst_name=""
            for x,y in systems.iteritems():
                if y>=max_syst:
                    tmp_syst_name=x
                    max_syst=y
        
            self.putative_system = tmp_syst_name
            self.systems=systems
            
            if len(systems.keys()) == 1:
                self.state = "clear"
            else:
                self.state = "ambiguous"


def build_clusters(hits, cfg):
    """
    Gets sets of contiguous hits according to the minimal inter_gene_max_space between two genes. Only for \"ordered\" datasets. 
    Remains to do : 
    
    - Implement case of circular replicons => need to store max position for each replicon ! and update... 
    - Runs on data grouped by replicon : does not check that genes from different replicons are not aggregated in a cluster (but functions Cluster and ClustersHandler do) 
    
    :param hits: a list of Hmmer hits to analyze 
    :type hits: a list of :class:`txsscanlib.report.Hit`
    :param cfg: the configuration object built from default and user parameters.
    :type cfg: :class:`txsscanlib.config.Config`
    :return: a set of clusters
    :rtype: :class:`txsscanlib.search_systems.ClustersHandler`
    """    
    
    _log.debug("Starting cluster detection with build_clusters... ")
    
    # Deals with different dataset types using Pipeline ?? 
    clusters = ClustersHandler(cfg)
    prev = hits[0]
    cur_cluster = Cluster(cfg)
    positions = []
    loner_state=False
    
    tmp=""
    for cur in hits[1:]:
        
        prev_max_dist = prev.get_syst_inter_gene_max_space()
        cur_max_dist = cur.get_syst_inter_gene_max_space()
        inter_gene = cur.get_position() - prev.get_position() - 1
        
        tmp="\n****\n"
        tmp+="prev_max_dist : %d\n"%(prev_max_dist)
        tmp+="cur_max_dist : %d\n"%(cur_max_dist)
        tmp+="Intergene space : %d\n"%(inter_gene)
        tmp+="Cur : %s"%cur
        tmp+="Prev : %s"%prev
        tmp+="Len cluster: %d\n"%len(cur_cluster)
        
        #print tmp
        
        # First condition removes duplicates (hits for the same sequence)
        # the two others takes into account either system1 parameter or system2 parameter    
        if(inter_gene <= prev_max_dist or inter_gene <= cur_max_dist ):
        # First check the cur.id is different from  the prev.id !!!
        #if(inter_gene!=-1 and (inter_gene <= prev_max_dist or inter_gene <= cur_max_dist )):
            #print "zero"
            if positions.count(prev.position) == 0:
                #print "un - ADD prev in cur_cluster"
                cur_cluster.add(prev)
                positions.append(prev.position)
                
            if positions.count(cur.position) == 0:
                #print "deux - ADD cur in cur_cluster"
                cur_cluster.add(cur)
                positions.append(cur.position)
                
            if prev.gene.loner:
                # PB !!!! Do not enter here when T3SS sctC loner and T2SS searched too... CHECK !!
                #print "trois - loner_state"
                #print "--- PREVLONER %s %s"%(prev.id, prev.gene.name)
                loner_state = True
        else:
            # Storage of the previous cluster
            if len(cur_cluster)>1:
                #print "quatre - ADD cur_cluster"
                clusters.add(cur_cluster) # Add an in-depth copy of the object? 
                #print(cur_cluster)
                
                cur_cluster = Cluster(cfg) 
                loner_state = False
                
            elif len(cur_cluster) == 1 and loner_state == True: # WTF?
                #print "cinq - ADD cur_cluster"
                #print "LNOER??"
                #print "PREVLONER %s %s"%(prev.id, prev.gene.name)
                clusters.add(cur_cluster) # Add an in-depth copy of the object? 
                #print(cur_cluster)
                    
                cur_cluster = Cluster(cfg) 
                loner_state = False
            
            if prev.gene.loner:
                #print "six - check"
                #print "PREVLONER ?? %s %s"%(prev.id, prev.gene.name)
                
                if positions.count(prev.position) == 0:
                    #print "six - ADD prev in cur_cluster, ADD cur_cluster"
                    cur_cluster.add(prev) 
                    clusters.add(cur_cluster) 
                    #print(cur_cluster)
                    
                    positions.append(prev.position)
                    loner_state = False  
                    cur_cluster = Cluster(cfg)  
                
            cur_cluster = Cluster(cfg)     
            
        prev=cur
    
    # !!! Treat the last cluster?     
    return clusters


class SystemOccurence(object):
    """
    This class is instantiated for a system that has been asked for detection. It can be filled step by step with hits. 
    A decision can then be made according to parameters defined *e.g.* quorum of genes.     
    """
    def __init__(self, system):
        """
        :param system: the system to \"fill\" with hits.
        :type system: :class:`txsscanlib.system.System` 
        """
        self.system = system
        self.system_name = system.name
        # Make those attributes non modifiable?
        self.mandatory_genes = {}
        for g in system.mandatory_genes:
            self.mandatory_genes[g.name]=0
        
        self.allowed_genes = {}
        for g in system.allowed_genes:
            self.allowed_genes[g.name]=0
        
        self.forbidden_genes = {}
        for g in system.forbidden_genes:
            self.forbidden_genes[g.name]=0
        
        #self.mandatory_genes = system.mandatory_genes
        #self.allowed_genes = system.allowed_genes
        #self.forbidden_genes = system.forbidden_genes
                 
    def __str__(self):
        out=""
        out+="Mandatory\n"
        for k, g in self.mandatory_genes.iteritems():
            out+="%s\t%d\n"%(k, g) 
        out+="Allowed\n"  
        for k, g in self.allowed_genes.iteritems():
            out+="%s\t%d\n"%(k, g)  
        out+="Forbidden\n"  
        for k, g in self.forbidden_genes.iteritems():
            out+="%s\t%d\n"%(k, g)   
        return out
        
    def add(self, cluster):
        
        for hit in cluster.hits:
            if hit.gene.is_mandatory(self.system):
                self.mandatory_genes[hit.gene.name]+=1
                #print "mandat"
            elif hit.gene.is_allowed(self.system):
                self.allowed_genes[hit.gene.name]+=1
                #print "allowed"
            elif hit.gene.is_forbidden(self.system):
                self.forbidden_genes[hit.gene.name]+=1
                #print "forbid"
            else:
                print "Foreign gene %s in cluster %s"%(hit.gene.name, self.system_name)
                
    # Add function in Genes isMandatory, isForbidden, isAllowed, isHomolog? 
    #def     

def analyze_clusters_replicon(clusters, systems, cfg):
    """
    Analyzes sets of contiguous hits (clusters) stored in a ClustersHandler for system detection:
        
    - split clusters if needed
    - delete them if they are not relevant
    - check the QUORUM for each system to detect, *i.e.* mandatory + allowed - forbidden
          
    Only for \"ordered\" datasets representing a whole replicon. 
    Reports systems occurence. 
    
    :param clusters: the set of clusters to analyze
    :type clusters: :class:`txsscanlib.search_systems.ClustersHandler` 
    :param systems: the set of systems to detect
    :type systems: a list of :class:`txsscanlib.system.System`
    :return: a set of systems occurence filled with hits found in clusters
    :rtype: a list of :class:`txsscanlib.search_systems.SystemOccurence` 
    
    """
    
    systems_occurences={}
    #for system in systems:
    #    systems_occurences[system.name]=SystemOccurence(system)
    
    syst_dict={}
    for system in systems:
        syst_dict[system.name]=system
    
    for clust in clusters.clusters:
        if clust.state == "clear":
            #so=systems_occurences[clust.putative_system]
            
            so=SystemOccurence(syst_dict[clust.putative_system])
            so.add(clust)
            print clust
            print so
        
        #break



def search_systems(hits, systems, cfg):
    """
    Runs systems search from hits according to the kind of database provided. 
    Analyze quorum and colocalization if required for system detection. 
    """

    if cfg.db_type == 'gembase':
        #build_clusters(hits, cfg)
        
        # Use of the groupby() function from itertools : allows to group Hits by replicon_name, 
        # and then apply the same build_clusters functions to replicons from "gembase" and "ordered_replicon" types of databases.
                
        #build_clusters(sub_hits, cfg) for sub_hits in [list(g) for k, g in itertools.groupby(hits, operator.attrgetter('replicon_name'))]
        for k, g in itertools.groupby(hits, operator.attrgetter('replicon_name')):
            sub_hits=list(g)
            print "\n************\nBuilding clusters for %s \n************\n"%k
            clusters=build_clusters(sub_hits, cfg)
            print "\n************\nAnalyzing clusters for %s \n************\n"%k
            analyze_clusters_replicon(clusters, systems, cfg)
            #break
            
    elif cfg.db_type == 'ordered_replicon':
        clusters=build_clusters(hits, cfg)
    elif cfg.db_type == 'unordered_replicon':
        pass
    elif cfg.db_type == 'unordered':
        pass
    else:
        raise ValueError("Invalid database type. ")





