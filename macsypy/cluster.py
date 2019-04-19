import itertools


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
        if 0 < dist <= model.max_gene_inter_space:
            return True
        elif dist < 0 and rep_info.topology == 'circular':
            # h1 et h2 sont par et d'autres ori
            dist = rep_info.length - h2.get_position + h1.get_position()
            return dist <= model.max_gene_inter_space
        return False


    clusters = []
    cluster_building = []
    # trier les hits par position croissantes
    hits.sort(key=lambda hit: hit.position_hit)
    # enlever les duplicats (meme sequence)
    hits = ([next(group) for pos, group in itertools.groupby(hits, lambda h: h.position)])

    cluster_building.append(hits[0])
    previous_hit = hits[0]
    for hit in hits:
        if is_in_cluster(previous_hit, hit):
            cluster_building.append(hit)
        else:
            if len(cluster_building) > 1:
                cluster = Cluster(cluster_building, model)
                clusters.append(cluster)
            cluster_building = []
        previous_hit = hit

    if len(cluster_building) > 1:
        # creer un cluster
        # si dernier hit du dernier cluster et 1er hit du 1er cluster < is_in_cluster
        #   fusionner les 2 clusters  dernier + 1er
        #sinon ajouter dernier cluster a la liste
    else:
        # regarder si previous hit et 1er hit du 1er cluster < is_in_cluster
        # ajouter hit au 1er cluster



    # ajouter le 1er elt de la liste dans une liste potentielle
    # calculer la distance entre 2 hits (en temr de nbre de genes)
    # si cette distance <= max_gene_inter_space
    #   ajouter le hit dans la liste du cluster
    # si non arreter la liste
    #    si la liste < 2 jeter la liste

    # cas des loner
    # cas des replicons circulaire
    #   => fusion du 1er et dernier cluster
    #  extension des derniers hits avec les premiers

    # replicon or RepliconInfo
    # Replicon info contains
    #  topology
    #  _min: rank of the first sequence in file
    #  _max: rank of the first sequence in file
    #  genes: (seq_name, seq_len)
    # taille du replicon ?
    #    en nombre de genes OK
    #    en nombre d'aa ??


class Cluster:

    def __init__(self, hits, model):
        self.hits = hits
        self.model = model


    def match(self, model):
        """
        test if this cluster match a model
        matching a model means satisfy the quorum constraints:
          * number of mandatory,
          * accessory,
          * presence of forbidden
        :param model: the model to chalenge
        :type model: :class:`macsypy.model.Model` object
        :return: a Putative System
        :rtype: :class:`Putative System` object
        """
        return True