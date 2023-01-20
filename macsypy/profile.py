#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2023  Institut Pasteur (Paris) and CNRS.           #
# See the COPYRIGHT file for details                                    #
#                                                                       #
# This file is part of MacSyFinder package.                             #
#                                                                       #
# MacSyFinder is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# MacSyFinder is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details .                         #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with MacSyFinder (COPYING).                                     #
# If not, see <https://www.gnu.org/licenses/>.                          #
#########################################################################

"""
Manage HMM profiles and hmmsearch execution
"""

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
                path = model_location.get_profile(gene.name)
            except KeyError:
                raise MacsypyError(f"'{model_location.name}/{gene.name}': No such profile")
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
        self.len, self.ga_threshold = self._profile_features()
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

    def _profile_features(self):
        """
        Parse the HMM profile to extract the length and the presence of GA bit threshold

        :return: the length, presence of ga bit threshold
        :rtype: tuple(int length, bool ga_threshold)
        """
        length = None
        ga_threshold = False
        with open(self.path) as hmm_file:
            for line in hmm_file:
                if line.startswith('LENG'):
                    length = int(line.split()[1])
                elif line.startswith('GA'):
                    try:
                        header, thld_1, thld_2 = line.split()
                    except ValueError:
                        _log.warning(f"{self.gene.name} GA score is not well formatted. expected: "
                                     f"'GA float float' got '{line.rstrip()}'.")
                        _log.warning(f"GA score will not used for gene '{self.gene.name}'.")
                        continue
                    if thld_2.endswith(';'):
                        thld_2 = thld_2[:-1]
                    try:
                        thld_1 = float(thld_1)
                        thld_2 = float(thld_2)
                        ga_threshold = True
                    except ValueError:
                        _log.warning(f"{self.gene.name} GA score is not well formatted "
                                     f"expected 2 floats got '{thld_1}' '{thld_2}'.")
                        _log.warning(f"GA score will not used for gene '{self.gene.name}'.")
                        continue
                elif line.startswith('STATS LOCAL'):
                    break
        return length, ga_threshold


    def __str__(self):
        """
        Print the name of the corresponding gene and the path to the HMM profile.
        """
        return f"{self.gene.name} : {self.path}"


    def execute(self, cpu=1):
        """
        Launch the Hmmer search (hmmsearch executable) with this profile

        :param int cpu: the number of cpu to use for hmmsearch (must be >= 1)
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
                if not self.cfg.cut_ga():
                    hmmer_threshold = f"-E {self.cfg.e_value_search():f}"
                elif self.cfg.cut_ga() and self.ga_threshold:
                    hmmer_threshold = "--cut_ga"
                else:
                    # cut_ga is True set but there is not self.ga_threshold:
                    hmmer_threshold = f"-E {self.cfg.e_value_search():f}"
                    _log.warning(f"GA bit thresholds unavailable on profile {self.gene.name}. "
                                 f"Switch to e-value threshold ({hmmer_threshold})")

                command = f'"{self.cfg.hmmer()}" --cpu {cpu} -o "{output_path}" {hmmer_threshold} ' \
                          f'"{self.path}" "{self.cfg.sequence_db()}" '
                _log.debug(f"{self.gene.name} Hmmer command line : {command}")
                try:
                    hmmer = Popen(command,
                                  shell=True,
                                  stdout=None,
                                  stdin=None,
                                  stderr=err_file,
                                  close_fds=False,
                                  )
                except Exception as err:
                    msg = f"Hmmer execution failed: command = {command} : {err}"
                    _log.critical(msg, exc_info=True)
                    raise err
                hmmer.wait()

            if hmmer.returncode != 0:
                if hmmer.returncode == -15:
                    msg = f"The Hmmer execution was aborted: command = {command} : " \
                          f"return code = {hmmer.returncode:d} check {err_path}"
                    _log.critical(msg)
                    return
                else:
                    msg = f"an error occurred during Hmmer execution: command = {command} : " \
                          f"return code = {hmmer.returncode:d} check {err_path}"
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
