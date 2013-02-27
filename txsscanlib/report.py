# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 29, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import os
import logging
_log = logging.getLogger('txsscan.' + __name__)


import abc
from threading import Lock

class HMMReport(object):
    """
    handle HMM report. extract a synthetic report from the raw hmmer output
    this class is an abstract class. there is 2 implementation of this abstract class
    depending if the genome baes is ordered or not.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, gene, hmmer_output, cfg):
        """
        :param gene: the gene corresponding to this profile
        :type gene: :class:`txsscanlib.gene.Gene` object
        :param hmmer_output: The path to hmmer output file
        :type hmmer_output: string
        :param cfg: the configuration 
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.gene = gene
        self._hmmer_raw_out = hmmer_output
        self._extract_out = None
        self.hits = []
        self.cfg = cfg
        self._lock = Lock()

    @abc.abstractmethod
    def extract(self):
        """
        Parse the output file of hmmer and produced a new synthetic report file.
        containing selected and sorted hits. ( **abstract method must be implemented in inherited classes** )
        """
        pass

    def __str__(self):
        s = "# gene: %s extract from %s hmm output\n" % (self.gene.name, self._hmmer_raw_out)
        s += "# profile length= %d\n" % len(self.gene.profile)
        s += "# i_evalue threshold= %f\n" % self.cfg.i_evalue_sel
        s += "# coverage threshold= %f\n" % self.cfg.coverage_profile
        #s += "# hit_id replicon_name position_hit gene_name gene_system i_eval score coverage\n"
        s += "# hit_id replicon_name position_hit gene_name gene_system i_eval score profile_coverage sequence_coverage begin end\n"
        for h in self.hits:
            s += str(h)
        return s

    def save_extract(self):
        """
        write the string representation of the extact report in a file.
        the name of this file is the concatenation of the gene and the res_extract_suffix from config
        """
        with self._lock:
            extract_out_name = self.gene.name + self.cfg.res_extract_suffix
            self._extract_out = os.path.join(self.cfg.working_dir, extract_out_name)
            with open(self._extract_out, 'w') as _file:
                _file.write(str(self))

    def best_hit(self):
        """
        return the best hit when multiple hits
        """
        try:
            return self.hits[0]
        except IndexError:
            return None

class UnOrderedHMMReport(HMMReport):
    """
    handle HMM report. extract a synthetic report from the raw hmmer output
    """

    pass


class OrderedHMMReport(HMMReport):
    """
    handle HMM report. extract a synthetic report from the raw hmmer output
    """


    def extract(self):
        """
        Parse the output file of hmmer compute from an unordered genes base 
        and produced a new synthetic report file.
        """
        with self._lock:
            # so the extract of a given HMM output is executed only once per run
            # if this method is called several times the first call induce the parsing of HMM out
            # the other calls do nothing
            if self.hits :
                return 
            i_evalue_sel = self.cfg.i_evalue_sel
            coverage_treshold = self.cfg.coverage_profile
            gene_profile_lg = len(self.gene.profile)
            with open(self._hmmer_raw_out, 'r') as hmm_out:
                for line in hmm_out:
                    if line.startswith(">> "):
                        fields = line.split()
                        hit_id = line.split()[1]
                        fields_hit = hit_id.split('_')
                        replicon_name = fields_hit[0]
                        position_hit = int(fields_hit[1]) / 10
                        # skip next 2 line
                        # the hits begins on the 3rd line
                        for _ in range(3):
                            line = hmm_out.next()
                        while not line.startswith("  Alignments"):
                            fields = line.split()
                            if(len(fields) > 1 and float(fields[5]) <= i_evalue_sel):
                                cov_profile = (float(fields[7]) - float(fields[6]) + 1) / gene_profile_lg
                                begin = int(fields[9])
                                end = int(fields[10])
                                cov_gene = (end - begin +1) #/ self.gene.sequence_lg # To be added in Gene: sequence_length
                                if (cov_profile >= coverage_treshold):
                                    i_eval = float(fields[5])
                                    score = float(fields[2])
                                    self.hits.append(Hit(self.gene,
                                                         self.gene.system,
                                                         hit_id,
                                                         replicon_name,
                                                         position_hit,
                                                         i_eval,
                                                         score,
                                                         cov_profile, 
                                                         cov_gene, 
                                                         begin, 
                                                         end))
                            line = hmm_out.next()
                self.hits.sort()


class Hit(object):
    """
    handle hits found by HMM. the hit are instanciate by :py:meth:`HMMReport.extract` method
    """
    
    #def __init__(self, gene, system, hit_id, replicon_name, position_hit, i_eval, score, coverage):
    def __init__(self, gene, system, hit_id, replicon_name, position_hit, i_eval, score, profile_coverage, sequence_coverage, begin_match, end_match):
        """
        :param gene: the gene corresponding to this profile
        :type gene: :class:`txsscanlib.gene.Gene` object
        :param system: the system to which this gene belongs
        :type system: :class:`txsscanlib.system.System` object
        :param hit_id: the identifier of the hit
        :type hit_id: string
        :param replicon_name: the name of the replicon
        :type replicon_name: string
        :param position_hit: the position of the hit on the sequence?
        :type position_hit: integer
        :param i_eval: the best-domain evalue
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
        self.system = system
        self.id = hit_id
        self.replicon_name = replicon_name
        self.position = position_hit
        self.i_eval = i_eval
        self.score = score
        self.profile_coverage = profile_coverage
        self.sequence_coverage = sequence_coverage
        self.begin_match = begin_match
        self.end_match = end_match

    def __str__(self):
        return "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t%d\n" % (self.id,
                                                     self.replicon_name,
                                                     self.position,
                                                     self.gene.name,
                                                     self.system.name,
                                                     self.i_eval,
                                                     self.score,
                                                     self.profile_coverage, 
                                                     self.sequence_coverage,
                                                     self.begin_match,
                                                     self.end_match)
                                                     
    def __cmp__(self, other):
        if self.id == other.id:
            if not self.gene.is_homolog(other.gene): 
                _log.warning("Non homologs match: %s (%s) %s (%s) for %s"%(self.gene.name, self.system.name, other.gene.name, other.system.name, self.id))
            return cmp(self.score, other.score)
        else:
            return cmp(self.id, other.id)
 
    def __eq__(self, other):
        return ( 
                self.gene.name == other.gene.name and
                self.system.name == other.system.name and
                self.id == other.id and
                self.replicon_name == other.replicon_name and
                self.position == other.position and
                self.i_eval == other.i_eval and
                self.score == other.score and
                self.profile_coverage == other.profile_coverage and
                self.sequence_coverage == other.sequence_coverage and
                self.begin_match == other.begin_match and
                self.end_match == other.end_match 
                )




 
