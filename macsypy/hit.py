import logging
_log = logging.getLogger(__name__)


class Hit(object):
    """
    Handle the hits filtered from the Hmmer search. The hits are instanciated by :py:meth:`HMMReport.extract` method
    """


    def __init__(self, gene, model, hit_id, hit_seq_length, replicon_name,
                 position_hit, i_eval, score, profile_coverage, sequence_coverage, begin_match, end_match):
        """
        :param gene: the gene corresponding to this profile
        :type gene: :class:`macsypy.gene.Gene` object
        :param model: the model to which this gene belongs
        :type model: :class:`macsypy.model.Model` object
        :param hit_id: the identifier of the hit
        :type hit_id: string
        :param hit_seq_length: the length of the hit sequence
        :type hit_seq_length: integer
        :param replicon_name: the name of the replicon
        :type replicon_name: string
        :param position_hit: the rank of the sequence matched in the input dataset file
        :type position_hit: integer
        :param i_eval: the best-domain evalue (i-evalue, "independent evalue")
        :type i_eval: float
        :param score: the score of the hit
        :type score: float
        :param profile_coverage: percentage of the profile that matches the hit sequence
        :type profile_coverage: float
        :param sequence_coverage: percentage of the hit sequence that matches the profile
        :type sequence_coverage: float
        :param begin_match: where the hit with the profile starts in the sequence
        :type begin_match: integer
        :param end_match: where the hit with the profile ends in the sequence
        :type end_match: integer
        """
        self.gene = gene
        self.model = model
        self.id = hit_id
        self.seq_length = hit_seq_length
        self.replicon_name = replicon_name
        self.position = position_hit
        self.i_eval = i_eval
        self.score = score
        self.profile_coverage = profile_coverage
        self.sequence_coverage = sequence_coverage
        self.begin_match = begin_match
        self.end_match = end_match


    def __str__(self):
        """
        Print useful information on the Hit: regarding Hmmer statistics, and sequence information
        """
        return "{id}\t{replicon_name}\t{position:d}\t{seq_len:d}\t{gene_name}\t{model_name}\t{i_evalue:.3e}\t{score:.3f}\
\t{profil_cov:.3f}\t{seq_cov:.3f}\t{begin_match:d}\t{end_match:d}\n".format(id=self.id,
                                                                            replicon_name=self.replicon_name,
                                                                            position=self.position,
                                                                            seq_len=self.seq_length,
                                                                            gene_name=self.gene.name,
                                                                            model_name=self.model.name,
                                                                            i_evalue=self.i_eval,
                                                                            score=self.score,
                                                                            profil_cov=self.profile_coverage,
                                                                            seq_cov=self.sequence_coverage,
                                                                            begin_match=self.begin_match,
                                                                            end_match=self.end_match)


    def __lt__(self, other):
        """
        compare two Hits. If the sequence identifier is the same, do the comparison on the score.
        Otherwise, do it on alphabetical comparison of the sequence identifier.

        :param other: the hit to compare to the current object
        :type other: :class:`macsypy.report.Hit` object
        :return: True if self is < other, False otherwise
        """
        if self.id == other.id:
            if not self.gene.is_homolog(other.gene):
                _log.warning("Non homologs match: {g_name} ({model_name}) {other_g_name} "
                             "({other_mod_name}) for {id}".format(g_name=self.gene.name,
                                                                  model_name=self.model.name,
                                                                  other_g_name=other.gene.name,
                                                                  other_mod_name=other.model.name,
                                                                  id=self.id))
            return self.score < other.score
        else:
            return self.id < other.id


    def __gt__(self, other):
        """
        compare two Hits. If the sequence identifier is the same, do the comparison on the score.
        Otherwise, do it on alphabetical comparison of the sequence identifier.

        :param other: the hit to compare to the current object
        :type other: :class:`macsypy.report.Hit` object
        :return: True if self is > other, False otherwise
        """
        if self.id == other.id:
            if not self.gene.is_homolog(other.gene):
                _log.warning("Non homologs match: {g_name} ({model_name}) {other_g_name} "
                             "({other_mod_name}) for {id}".format(g_name=self.gene.name,
                                                                  model_name=self.model.name,
                                                                  other_g_name=other.gene.name,
                                                                  other_mod_name=other.model.name,
                                                                  id=self.id))
            return self.score > other.score
        else:
            return self.id > other.id


    def __eq__(self, other):
        """
        Return True if two hits are totally equivalent, False otherwise.

        :param other: the hit to compare to the current object
        :type other: :class:`macsypy.report.Hit` object
        :return: the result of the comparison
        :rtype: boolean
        """
        epsilon = 0.001
        return (self.gene.name == other.gene.name and
                self.model.name == other.model.name and
                self.id == other.id and
                self.seq_length == other.seq_length and
                self.replicon_name == other.replicon_name and
                self.position == other.position and
                abs(self.i_eval - other.i_eval) <= epsilon and
                abs(self.score - other.score) <= epsilon and
                abs(self.profile_coverage - other.profile_coverage) <= epsilon and
                abs(self.sequence_coverage - other.sequence_coverage) <= epsilon and
                self.begin_match == other.begin_match and
                self.end_match == other.end_match
                )


    def get_position(self):
        """
        :returns: the position of the hit (rank in the input dataset file)
        :rtype: integer
        """
        return self.position


    def get_syst_inter_gene_max_space(self):
        """
        :returns: the 'inter_gene_max_space' parameter defined for the gene of the hit
        :rtype: integer
        """
        return self.gene.model.inter_gene_max_space


class HitRegistry:

    def __init__(self):
        self._DB = {}

    def __contains__(self, hit):
        return hit in self._DB

    def __getitem__(self, hit):
        return self._DB[hit]

    def __setitem__(self, hit, system):
        if hit in self._DB:
            self._DB[hit].append(system)
        else:
            self._DB[hit] = [system]

