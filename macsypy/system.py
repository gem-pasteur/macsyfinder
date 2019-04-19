import abc
from itertools import chain
from collections import Counter


def match(cluster, model):
    def create_exchangeable_map(genes):
        map = {}
        for gene in model.mandatory_genes:
            for ex_gene in chain(gene.get_homolgs(), gene.get_analogs()):
                map[ex_gene.name] = gene.name
        return map

    mandatory_counter = {g.name: 0 for g in model.mandatory_genes}
    exchangeable_mandatory = create_exchangeable_map(model.mandatory_genes)

    accessory_counter = {g.name: 0 for g in model.accessory_genes}
    exchangeable_accessory = create_exchangeable_map(model.forbidden_genes)

    forbidden_counter = {g.name: 0 for g in model.forbidden_genes}
    exchangeable_forbidden = create_exchangeable_map(model.forbidden_genes)

    # Forbidden genes do not play a role in the system, thus they do not have the multi_system feature
    multi_syst_genes = {g.name: 0 for g in chain(model.mandatory_genes, model.forbidden_genes) if g.multi_system}

    genes_occ = Counter([h.gene for h in cluster.hits])

    for gene, occ in genes_occ.items():
        if gene.name in mandatory_counter :
            mandatory_counter[gene.name] += occ
        elif gene.name in exchangeable_mandatory:
            gene_ref = exchangeable_mandatory[gene.name]
            mandatory_counter[gene_ref] += occ
        elif gene.name in accessory_counter or gene.name in :
            accessory_counter[gene.name] += occ
        elif gene.name in exchangeable_accessory:
            gene_ref = exchangeable_accessory[gene.name]
            accessory_counter[gene_ref] += occ
        elif gene.name in forbidden_counter:
            forbidden_counter[gene.name] += occ
        elif gene.name in exchangeable_forbidden:
            gene_ref = exchangeable_forbidden[gene.name]
            forbidden_counter[gene_ref] += occ

    # a la fin
    # nombre de gene_mandatory > min_mandatory_gene_required?
    # min_mandatory_gene_required en en nombre d'occurence ou nombre de genes differents
    # nombre de gene accessory > min_gene_required ?
    # min_gene_required c'est accessory ou accessory + mandatory ?
    # nombre de gene forbiden
    # si satifait contraintes, instancie un putative system
    # sinon cree un RejectedCluster



class System(metaclass=abc.ABCMeta):

    def __init__(self, model, hits, mandatory_counter, accessory_counter, forbiden_counter,
                 exchangeable_mandatory, exchangeable_accessory, exchangeable_forbiden, multi_syst_genes):
        self.model = model
        self.hits = hits
        self.mandatory_counter
        self.accessory_counter
        self.forbiden_counter



    @abc.abstractmethod
    def __str__(self):
        """

        :return:
        """
        return ""

    @abc.abstractmethod
    def to_json(self):
        """

        :return:
        """
        return {}


class PutativeSystem(System):

    def __str__(self):
        s = ""
        if self.mandatory_genes:
            for g_name, g_occ in self.mandatory_genes.items():
                s += "{}\t{:d}\n".format(g_name, g_occ)
        if self.accessory_genes:
            for g_name, g_occ in self.accessory_genes.items():
                s += "{}\t{:d}\n".format(g_name, g_occ)
        if self.forbidden_genes:
            for g_name, g_occ in self.forbidden_genes.items(g_name, g_occ):
                s += "{0}\t{1:d}\n".format(k, g)

        if self.multi_syst_genes:
            for g_name, g_occ in self.multi_syst_genes.items():
                s += "{0}\t{:d}\n".format(g_name, g_occ)
        return s


class ConfirmedSystem(System):
    pass

class AmbiguousSystem(System):
    pass

class RejectedSystem(System):
    pass