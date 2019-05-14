import itertools
import json

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
    print("################################################# DEBUT MATCH ###########################################")
    print("########### system match L36 clusters", clusters)
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
    print("############# ENTER LOOP on CLUSTER ################################")
    for cluster in clusters:
        print("========================")
        print("cluster", cluster)
        valid_hits = []
        for hit in cluster.hits:
            print("----------------------)")
            print("hit", hit)
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
        print("######### =============================")
        print("mandatory_counter", mandatory_counter)
        print("accessory_counter", accessory_counter)
        print("forbidden_counter", forbidden_counter)
        print("######### =============================")
        if valid_hits:
            valid_clusters.append(Cluster(valid_hits, model))
    # the count is finished
    # check if the quorum is reached
    # count how many different genes are represented in the clusters
    mandatory_genes = [g for g, occ in mandatory_counter.items() if occ > 0]
    accessory_genes = [g for g, occ in accessory_counter.items() if occ > 0]
    forbidden_genes = [g for g, occ in forbidden_counter.items() if occ > 0]
    print("######### =============================")
    print("mandatory_genes", mandatory_genes)
    print("accessory_genes", accessory_genes)
    print("forbidden_genes", forbidden_genes)
    print("######### =============================")
    reasons = []
    is_a_system = True
    if forbidden_genes:
        is_a_system = False
        reasons.append('There is {} forbidden genes occurrence(s): {}'.format(
            len(forbidden_hits), ', '.join(h.gene.name for h in forbidden_hits)
        ))
    if len(mandatory_genes) < model.min_mandatory_genes_required:
        is_a_system = False
        reasons.append('The quorum of mandatory genes required ({}) is not reached: {}'.format(
            model.min_mandatory_genes_required, len(mandatory_genes)))
    if len(accessory_genes) + len(mandatory_genes) < model.min_genes_required:
        is_a_system = False
        reasons.append('The quorum of genes required ({}) is not reached: {}'.format(
            model.min_genes_required, len(accessory_genes)
        ))

    if is_a_system:
        res = System(model, valid_clusters)
        for hit in valid_hits:
            if hit.gene.multi_system:  # gene or gene_ref ?
                hit_registry[hit] = res
    else:
        reason = '\n'.join(reasons)
        res = RejectedClusters(model, clusters, reason)
    print("############## system match L130 res \n", res, type(res))
    print("##################################### FIN system match ########################################################")
    return res, hit_registry


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
        self._mandatory_occ = {g.name: [] for g in self.model.mandatory_genes}
        self._accessory_occ = {g.name: [] for g in self.model.accessory_genes}
        for hit in self.hits:
            if hit.status == GeneStatus.MANDATORY:
                self._mandatory_occ[hit.gene_ref.name].append(hit)
            elif hit.status == GeneStatus.ACCESSORY:
                self._accessory_occ[hit.gene_ref.name].append(hit)

    @property
    def hits(self):
        hits = [h for cluster in self.clusters for h in cluster.hits]
        return hits


    @property
    def multi_loci(self):
        return len(self.clusters) > 1


    def __str__(self):

        s = """system id = {sys_id}
model = {model} 
loci nb = {loci}
replicon = {rep_name}
clusters = {clst}
""".format(sys_id=self.id,
           model=self.model.fqn,
           loci=len(self.clusters),
           rep_name=self._replicon_name,
           clst=", ".join(["[" + ", ".join([v_h.gene.name for v_h in cluster.hits]) + "]" for cluster in self.clusters])
           )
        for title, genes in (("mandatory", self._mandatory_occ), ("accessory", self._accessory_occ)):
            s += "\n{} genes:\n".format(title)
            for g_name, hits in genes.items():
                s += "\t- {g_ref}: {occ} ({hits})\n".format(g_ref=g_name,
                                                            occ=len(hits),
                                                            hits=', '.join([h.gene.name for h in hits])
                                                            )
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
