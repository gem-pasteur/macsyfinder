# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################

import itertools
import json
import statistics
from itertools import chain
import logging
_log = logging.getLogger(__name__)

from.gene import GeneStatus
from .cluster import Cluster, RejectedClusters
from .hit import ValidHit


# la liste des clusters a fournir est a generer avant match
# si len(clusters) = 1 single_loci
# si len(clusters) > 1 multi_loci
# il faut genegerer la liste de toutes les combinaisons
# et appeler cette fonction pour chaqu'une entre elles
# from itertools import combinations

# combinations('ABCD', 1) => inutile mais generique => single_loucs
# combinations('ABCD', 2) => multi_locus a ne faire que si model.multi_locus= True
# combinations('ABCD', 3)
# combinations('ABCD', len("ABCD")) => inutile mais generique => recheche parmis tous les clusters


def match(clusters, model, hit_registry):
    """
    Check a set of clusters fill model constraints.
    If yes create a :class:`macsypy.system.PutativeSystem` otherwise create
    a :class:`macsypy.cluster.RejectedClusters`.

    :param clusters: The list of cluster to check if fit the model
    :type clusters: list of :class:`macsypy.cluster.Cluster` objects
    :param model:  The model to consider
    :type model: :class:`macsypy.model.Model` object
    :param hit_registry: The registry where all hits => System are registered
    :return: either a System or a RejectedClusters
    :rtype: :class:`macsypy.system.System` or :class:`macsypy.cluster.RejectedClusters` object
    """
    def create_exchangeable_map(genes):
        """
        create a map between an exchangeable (homolog or analog) gene name and it's gene ref

        :param genes: The genes to get the homologs or analogs
        :type genes: list of :class:`macsypy.gene.Gene` objects
        :rtype: a dict with keys are the homolog_or analog gene_name the reference gene name
        """
        map = {}
        for gene in genes:
            if gene.exchangeable:
                for ex_gene in itertools.chain(gene.get_homologs(), gene.get_analogs()):
                    map[ex_gene.name] = gene
        return map

    # init my structures to count gene occurrences
    mandatory_counter = {g.name: 0 for g in model.mandatory_genes}
    exchangeable_mandatory = create_exchangeable_map(model.mandatory_genes)

    accessory_counter = {g.name: 0 for g in model.accessory_genes}
    exchangeable_accessory = create_exchangeable_map(model.accessory_genes)

    forbidden_counter = {g.name: 0 for g in model.forbidden_genes}
    exchangeable_forbidden = create_exchangeable_map(model.forbidden_genes)

    # count the hits
    # and track for each hit for which gene it counts for
    valid_clusters = []
    forbidden_hits = []
    for cluster in clusters:
        valid_hits = []
        for hit in cluster.hits:
            gene_name = hit.gene.name
            if gene_name in mandatory_counter:
                mandatory_counter[hit.gene.name] += 1
                valid_hits.append(ValidHit(hit, hit.gene, GeneStatus.MANDATORY))
            elif gene_name in exchangeable_mandatory:
                gene_ref = exchangeable_mandatory[gene_name]
                mandatory_counter[gene_ref.name] += 1
                valid_hits.append(ValidHit(hit, gene_ref, GeneStatus.MANDATORY))
            elif gene_name in accessory_counter:
                accessory_counter[gene_name] += 1
                valid_hits.append(ValidHit(hit, hit.gene, GeneStatus.ACCESSORY))
            elif gene_name in exchangeable_accessory:
                gene_ref = exchangeable_accessory[gene_name]
                accessory_counter[gene_ref.name] += 1
                valid_hits.append(ValidHit(hit, gene_ref, GeneStatus.ACCESSORY))
            elif gene_name in forbidden_counter:
                forbidden_counter[gene_name] += 1
                # valid_hits.append(ValidHit(hit, hit.gene, GeneStatus.FORBIDDEN))
                forbidden_hits.append(hit)
            elif gene_name in exchangeable_forbidden:
                gene_ref = exchangeable_forbidden[gene_name]
                forbidden_counter[gene_ref.name] += 1
                # valid_hits.append(ValidHit(hit, hit.gene.ref, GeneStatus.FORBIDDEN))
                forbidden_hits.append(hit)
        if valid_hits:
            valid_clusters.append(Cluster(valid_hits, model))
    # the count is finished
    # check if the quorum is reached
    # count how many different genes are represented in the clusters
    mandatory_genes = [g for g, occ in mandatory_counter.items() if occ > 0]
    accessory_genes = [g for g, occ in accessory_counter.items() if occ > 0]
    forbidden_genes = [g for g, occ in forbidden_counter.items() if occ > 0]
    _log.debug("#" * 50)
    _log.debug("mandatory_genes: {}".format(mandatory_genes))
    _log.debug("accessory_genes: {}".format(accessory_genes))
    _log.debug("forbidden_genes: {}".format(forbidden_genes))

    reasons = []
    is_a_system = True
    if forbidden_genes:
        is_a_system = False
        msg = 'There is {} forbidden genes occurrence(s): {}'.format(
            len(forbidden_hits), ', '.join(h.gene.name for h in forbidden_hits)
        )
        reasons.append(msg)
        _log.debug(msg)
    if len(mandatory_genes) < model.min_mandatory_genes_required:
        is_a_system = False
        msg = 'The quorum of mandatory genes required ({}) is not reached: {}'.format(
            model.min_mandatory_genes_required, len(mandatory_genes))
        reasons.append(msg)
        _log.debug(msg)
    if len(accessory_genes) + len(mandatory_genes) < model.min_genes_required:
        is_a_system = False
        msg = 'The quorum of genes required ({}) is not reached: {}'.format(
            model.min_genes_required, len(accessory_genes) + len(mandatory_genes)
        )
        reasons.append(msg)
        _log.debug(msg)

    if is_a_system:
        res = System(model, valid_clusters)
        for hit in valid_hits:
            if hit.gene.multi_system:  # gene or gene_ref ?
                hit_registry[hit] = res
        _log.debug("is a putative system")
    else:
        reason = '\n'.join(reasons)
        res = RejectedClusters(model, clusters, reason)
    _log.debug("#" * 50)
    return res, hit_registry


def track_multi_systems(systems):
    multi_sys_tracker = {}
    model_2_system = {}

    for system in systems:
        v_hits = system.hits

        model_fqn = system.model.fqn
        if model_fqn not in model_2_system:
            model_2_system[model_fqn] = set()
        model_2_system[model_fqn].add(system)

        for v_hit in v_hits:
            hit = v_hit.hit

            if hit not in multi_sys_tracker:
                multi_sys_tracker[hit] = set()
            multi_sys_tracker[hit].add(v_hit.gene_ref.model.fqn)

    for hit, models_fqn in multi_sys_tracker.items():
        for model_fqn in models_fqn:
            for system in model_2_system[model_fqn]:
                hit.add_system(system)


class System:

    _id = itertools.count(1)

    def __init__(self, model, clusters):
        """

        :param model:  The model which has ben used to build this system
        :type model: :class:`macsypy.model.Model` object
        :param clusters: The list of cluster that form this system
        :type clusters: list of :class:`macsypy.cluster.Cluster` objects
        """
        self._replicon_name = clusters[0].replicon_name
        self.id = "{}_{}_{}".format(self._replicon_name, model.name, next(self._id))
        self.model = model
        self.clusters = clusters
        self._mandatory_occ = None
        self._accessory_occ = None
        self._count()

    def _count(self):
        """
        fill 2 structures one for mandatory the other for accessory
        each structure count how many hit for each gene
        :return: None
        """
        self._mandatory_occ = {g.name: [] for g in self.model.mandatory_genes}
        self._accessory_occ = {g.name: [] for g in self.model.accessory_genes}

        # all the hits are ValidHit
        for hit in self.hits:
            if hit.status == GeneStatus.MANDATORY:
                self._mandatory_occ[hit.gene_ref.name].append(hit)
            elif hit.status == GeneStatus.ACCESSORY:
                self._accessory_occ[hit.gene_ref.name].append(hit)

    @property
    def wholeness(self):
        """

        :return:
        """
        # model completude
        score = sum([1 for hits in chain(self._mandatory_occ.values(), self._accessory_occ.values()) if hits]) / \
                (len(self._mandatory_occ) + len(self._accessory_occ))

        return score

    @property
    def score(self):
        """
        :return: a score take in account
            * if a hit match for the gene or is an homolog or analog
            * if a hit is duplicated and already present in the system or the cluster
            * if a hit match for mandatory/accessory gene of the model
        :rtype: float
        """
        weights = {'analog': 0.75,
                   'homolog': 0.75,
                   'in_clst': 0,
                   'in_syst': -0.5,
                   'mandatory': 1,
                   'accessory': 0.5,
                   }
        score = 0
        seen_in_system = set()
        clusters = sorted(self.clusters, key=len, reverse=True)
        for clst in clusters:
            seen_in_clst = set()
            for v_hit in clst.hits:
                # attribute a score for this hit according to mandatory/accessory
                if v_hit.gene == v_hit.gene_ref:
                    if v_hit.status == GeneStatus.MANDATORY:
                        hit_score = weights['mandatory']
                    elif v_hit.status == GeneStatus.ACCESSORY:
                        hit_score = weights['accessory']
                # weighted the hit score according to the hit match the gene or is an analog/homolog
                if v_hit.gene == v_hit.gene_ref:
                    pass
                elif v_hit.gene_ref.is_analog(v_hit.gene):
                    hit_score *= weights['analog']
                elif v_hit.gene_ref.is_homolog(v_hit.gene):
                    hit_score *= weights['homolog']

                # compute the global score by weighting the hit_score according if this gene is already present
                # in this cluster or the system
                if v_hit.gene_ref in seen_in_system:
                    hit_score *= weights['in_syst']
                elif v_hit.gene_ref in seen_in_clst:
                    hit_score *= weights['in_clst']
                else:
                    pass
                score += hit_score
                seen_in_clst.add(v_hit.gene_ref)
            seen_in_system.update(seen_in_clst)
        return score


    def occurence(self):
        """
        sometimes several systems collocates so they form only one cluster
        so macsyfinder build only one system
        the occurrence is an indicator of how many systems are
        it's based on the number of occurrence of each mandatory genes

        :return: a predict number of biologic systems
        """
        occ_per_gene = [len(hits) for hits in self._mandatory_occ.values()]
        # if a systems contains 5 gene whit occ of 1 and 5 gene with 0 occ
        # the median is 0.5
        # round(0.5) = 0
        # so I fix a floor value at 1
        return max(1, round(statistics.median(occ_per_gene)))


    @property
    def hits(self):
        """
        :return: The list of all hits that compose this system
        :rtype: [:class:`macsypy.hit.ValidHits` , ... ]
        """
        hits = [h for cluster in self.clusters for h in cluster.hits]
        return hits

    @property
    def loci(self):
        """
        :return: The number of loci of this system
        :rtype: int > 0
        """
        # we do not take loners in account
        loci = sum([1 for c in self.clusters if len(c) > 1])
        return loci


    @property
    def multi_loci(self):
        """
        :return: True if the systems is multi_loci. False otherwise
        :rtype: bool
        """
        return self.loci > 1


    def __str__(self):

        s = """system id = {sys_id}
model = {model}
replicon = {rep_name}
clusters = {clst}
occ = {occ}
wholeness = {wholeness:.3f}
loci nb = {loci}
score = {score:.3f}
""".format(sys_id=self.id,
           model=self.model.fqn,
           loci=self.loci,
           rep_name=self._replicon_name,
           clst=", ".join(["[" + ", ".join([str((v_h.gene.name, v_h.position)) for v_h in cluster.hits]) + "]"
                                                                               for cluster in self.clusters]),
           occ=self.occurence(),
           wholeness=self.wholeness,
           score=self.score
           )
        for title, genes in (("mandatory", self._mandatory_occ), ("accessory", self._accessory_occ)):
            s += "\n{} genes:\n".format(title)
            for g_name, hits in genes.items():
                s += "\t- {g_ref}: {occ} ".format(g_ref=g_name,
                                                  occ=len(hits))
                all_hits_str = []
                for h in hits:
                    used_in_systems = [s.id for s in h.used_in_systems() if s.model.fqn != self.model.fqn]
                    if used_in_systems:
                        hit_str = "{} [{}]".format(h.gene.name, ', '.join(used_in_systems))
                    else:
                        hit_str = "{}".format(h.gene.name)
                    all_hits_str.append(hit_str)
                s += "({})\n".format(", ".join(all_hits_str))

        return s


    def to_json(self):
        """
        :return: a serialisation of this system in json format
                 The json have the following structure
                 {'id': str system_id
                  'model': str model fully qualified name
                  'loci_nb': int number of loci
                  'replicon_name': str the replicon name
                  'clusters': [[ str hit gene name, ...], [...]]
                  'gene_composition': {
                        'mandatory': {str gene_ref name: [ str hit gene name, ... ]},
                        'accessory': {str gene_ref name: [ str hit gene name, ... ]}
                        }
                 }
        """
        system = {'id': self.id,
                  'model': self.model.fqn,
                  'loci_nb': len(self.clusters),
                  'replicon_name': self._replicon_name,
                  'clusters': [[v_h.gene.name for v_h in cluster.hits]for cluster in self.clusters],
                  'gene_composition':
                      {'mandatory': {gene_ref: [hit.gene.name for hit in hits]
                                     for gene_ref, hits in self._mandatory_occ.items()},
                       'accessory': {gene_ref: [hit.gene.name for hit in hits]
                                     for gene_ref, hits in self._accessory_occ.items()}
                       }
                  }
        return json.dumps(system)
