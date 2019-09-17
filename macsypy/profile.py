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


import os
import logging
_log = logging.getLogger(__name__)

from subprocess import Popen
from threading import Lock

from .report import GembaseHMMReport, GeneralHMMReport, OrderedHMMReport
from .error import MacsypyError


class ProfileFactory:
    """
    Build and store all Profile objects. Profiles must not be instanciated directly.
    The profile_factory must be used. The profile_factory ensures there is only one instance
    of profile for a given name.
    To get a profile, use the method get_profile. If the profile is already cached, this instance is returned.
    Otherwise a new profile is built, stored in the profile_factory and then returned.

    """

    def __init__(self, cfg):
        self._profiles = {}
        self.cfg = cfg

    def get_profile(self, gene, model_location):
        """
        :param gene: the gene associated to this profile
        :type gene: :class:`macsypy.gene.Gene` or :class:`macsypy.gene.Homolog` or :class:`macsypy.gene.Analog` object
        :param model_location: The where to get the profile
        :type model_location: :class:`macsypy.registries.ModelLocation` object.
        :return: the profile corresponding to the name.
                 If the profile already exists, return it. Otherwise build it, store it and return it.
        :rtype: :class:`macsypy.profile.Profile` object
        """
        key = (model_location.name, gene.name)
        if key in self._profiles:
            profile = self._profiles[key]
        else:
            try:
                path = model_location.get_profile(gene.name, )
            except KeyError:
                raise MacsypyError("{0}: No such profile".format(gene.name))
            profile = Profile(gene, self.cfg, path)
            self._profiles[key] = profile
        return profile


class Profile:
    """
    Handle a HMM protein profile
    """

    def __init__(self, gene, cfg, path):
        """

        :param gene: the gene corresponding to this profile
        :type gene: :class:`macsypy.secretion.Gene` object
        :param cfg: the configuration 
        :type cfg: :class:`macsypy.config.Config` object
        :param path: the path to the hmm profile.
        :type path: string
        """
        self.gene = gene
        self.path = path
        self.len = self._len()
        self.cfg = cfg 
        self.hmm_raw_output = None
        self._report = None
        self._lock = Lock()


    def __len__(self):
        """
        :return: the length of the HMM protein profile
        :rtype: int
        """
        return self.len


    def _len(self):
        """
        Parse the HMM profile file to get and store the length.
        This private method is called at the Profile init.
        """
        length = None
        with open(self.path) as f:
            for l in f:
                if l.startswith("LENG"):
                    length = int(l.split()[1])
                    break
        return length


    def __str__(self):
        """
        Print the name of the corresponding gene and the path to the HMM profile.
        """
        return "{0} : {1}".format(self.gene.name, self.path)


    def execute(self):
        """
        Launch the Hmmer search (hmmsearch executable) with this profile

        :return: an object storing information on the results of the HMM search (HMMReport)
        :rtype:  :class:`macsypy.report.HMMReport` object
        """
        with self._lock:
            # the results of HMM is cached 
            # so HMMsearch is executed only once per run
            # if this method is called several times,
            # the first call induce the execution of HMMsearch and generate a report
            # the other calls return directly this report
            if self._report is not None:
                return self._report
            hmmer_dir = os.path.join(self.cfg.working_dir(), self.cfg.hmmer_dir())
            if not os.path.exists(hmmer_dir):
                os.mkdir(hmmer_dir)
            output_path = os.path.join(hmmer_dir,  self.gene.name + self.cfg.res_search_suffix())
            err_path = os.path.join(hmmer_dir,
                                    self.gene.name + os.path.splitext(self.cfg.res_search_suffix())[0] + ".err")

            with open(err_path, 'w') as err_file:
                options = {"hmmer_exe": self.cfg.hmmer(),
                           "output_file": output_path,
                           "e_value_res": self.cfg.e_value_search(),
                           "profile": self.path,
                           "sequence_db": self.cfg.sequence_db(),
                           }
                command = "{hmmer_exe} --cpu 0 -o {output_file} -E {e_value_res:f} {profile} {sequence_db} "\
                    .format(**options)
                _log.debug("{0} Hmmer command line : {1}".format(self.gene.name, command))
                try:
                    hmmer = Popen(command,
                                  shell=True,
                                  stdout=None,
                                  stdin=None,
                                  stderr=err_file,
                                  close_fds=False,
                                  )
                except Exception as err:
                    msg = "Hmmer execution failed: command = {0} : {1}".format(command, err)
                    _log.critical(msg, exc_info=True)
                    raise err
                hmmer.wait()

            if hmmer.returncode != 0:
                if hmmer.returncode == -15:
                    msg = "The Hmmer execution was aborted: command = {0} : return code = {1:d} check {2}".format(
                        command, hmmer.returncode, err_path)
                    _log.critical(msg)
                    return
                else:
                    msg = "an error occurred during Hmmer execution: command = {0} : return code = {1:d} check {2}"\
                        .format(command, hmmer.returncode, err_path)
                    _log.debug(msg, exc_info=True)
                    _log.critical(msg)
                    raise RuntimeError(msg)
            self.hmm_raw_output = output_path
            db_type = self.cfg.db_type()
            if db_type == 'gembase':
                report = GembaseHMMReport(self.gene, output_path, self.cfg)
            elif db_type == 'ordered_replicon':
                report = OrderedHMMReport(self.gene, output_path, self.cfg)
            else:
                report = GeneralHMMReport(self.gene, output_path, self.cfg)
            self._report = report
            return report
