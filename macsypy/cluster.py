import itertools
import logging

from .error import MacsypyError

_log = logging.getLogger(__name__)


def build_clusters(hits, rep_info, model):
    """
    From a list of filtered hits, and replicon information (topology, length),
    build all lists of hits that satisfied the constraints:
        * max_gene_inter_space
        * lonner
        * multi_system
    if Yes create a cluster

    :param hits: list of filtered hits
    :type hits: list of :class:`macsypy.report.Hit` objects
    :param rep_info: the replicon to analyse
    :type rep_info: :class:`macsypy.Indexes.RepliconInfo` object
    :param model:
    :type model: :class:`macsypy.model.Model` object
    :return: list of clusters
    :rtype: List of :class:`Cluster` objects
    """
    def is_in_cluster(h1, h2):
        dist = h2.get_position() - h1.get_position() - 1
        inter_gene_max_space = max(h1.gene.inter_gene_max_space, h2.gene.inter_gene_max_space)
        if 0 < dist <= inter_gene_max_space:
            return True
        elif dist < 0 and rep_info.topology == 'circular':
            # h1 and h2 overlap the ori
            dist = rep_info.max - h1.get_position() + h2.get_position() - rep_info.min
            return dist <= inter_gene_max_space
        return False

    clusters = []
    cluster_scaffold = []
    # sort hits by increasing position and then descending score
    hits.sort(key=lambda h: (h.position, - h.score))
    # remove duplicates hits (several hits for the same sequence),
    # keep the first one, this with the best score
    # position == sequence rank in replicon
    hits = ([next(group) for pos, group in itertools.groupby(hits, lambda h: h.position)])
    cluster_scaffold.append(hits[0])
    previous_hit = hits[0]
    for hit in hits[1:]:
        if is_in_cluster(previous_hit, hit):
            cluster_scaffold.append(hit)
        else:
            if len(cluster_scaffold) > 1:
                cluster = Cluster(cluster_scaffold, model)
                clusters.append(cluster)
            cluster_scaffold = [hit]
        previous_hit = hit

    if len(cluster_scaffold) > 1:
        new_cluster = Cluster(cluster_scaffold, model)
        if is_in_cluster(new_cluster.hits[-1], clusters[0].hits[0]):
            clusters[0].merge(new_cluster, before=True)
        else:
            clusters.append(new_cluster)
    else:
        if is_in_cluster(previous_hit, clusters[0].hits[0]):
            clusters[0].merge(Cluster([previous_hit], model), before=True)

    return clusters


class Cluster:

    def __init__(self, hits, model):
        self.hits = hits
        self.model = model

    def merge(self, cluster, before=False):
        """
        merge the cluster in this one.

        :param cluster:
        :type cluster:
        :param bool before: If True the hits of the cluster will be add at the end of this one,
                            Otherwise the cluster hits be insert before the hits of this one
        :return: None
        :raise MasypyError: if the two cluster have not the same model
        """
        if cluster.model != self.model:
            raise MacsypyError("Try to merge Clusters from different model")
        else:
            if before:
                self.hits = cluster.hits + self.hits
            else:
                self.hits.extend(cluster.hits)


    def __str__(self):
        s = """Cluster:
    - model: {}
    - hits: {}""".format(self.model.name, ', '.join(["({}, {})".format(h.id,
                                                                       h.gene.name,
                                                                       h.position) for h in self.hits]))
        return s


class RejectedCluster:

    def __init__(self, cluster, reason):
        self.cluster = cluster
        self.reason = reason