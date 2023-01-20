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
Manage the hmm step (hmmsearch or recover results from previous run) in parallele
"""

import concurrent.futures
import logging
import signal
import sys
import shutil
import os.path
import math
from .report import GembaseHMMReport, GeneralHMMReport, OrderedHMMReport
from .utils import threads_available
from .database import Indexes

_log = logging.getLogger(__name__)


def worker_cpu(genes_nb, cfg):
    """
    Compute the optimum number of worker and cpu per worker
    The number of worker is set by the user (1 by default 0 means all worker available)

    we use one worker per gene
    if number of workers is greater than number of genes then several cpu can be use by
    hmsearch to speed up the search step

    :param int genes_nb: the number of genes to search
    :param cfg: The macsyfinder configuration
    :type cfg: :class:`macsypy.config.Config` object
    :return: the number of worker and cpu_per_worker to use
    :rtype: tuple (int worker_nb, int cpu_per_worker)
    """
    worker_nb = cfg.worker()
    if not worker_nb:
        worker_nb = threads_available()
    cpu_per_worker = max(1, math.floor(worker_nb / genes_nb))
    return worker_nb, cpu_per_worker


def search_genes(genes, cfg):
    """
    For each gene of the list, use the corresponding profile to perform an Hmmer search, and parse the output
    to generate a HMMReport that is saved in a file after CoreHit filtering.
    These tasks are performed in parallel using threads.
    The number of workers can be limited by worker_nb directive in the config object or
    in the command-line with the "-w" option.

    :param genes: the genes to search in the input sequence dataset
    :type genes: list of :class:`macsypy.gene.ModelGene` objects
    :param cfg: the configuration object
    :type cfg: :class:`macsypy.config.Config` object
    """
    worker_nb, cpu_per_worker = worker_cpu(len(genes), cfg)
    _log.debug(f"worker_nb = {worker_nb:d}\tcpu per worker = {cpu_per_worker}")
    all_reports = []

    def stop(signum, frame):
        """stop the main process, its threads and subprocesses"""
        _log.critical('KILL all Processes')
        proc_grp_id = os.getpgid(0)
        os.killpg(proc_grp_id, signum)
        sys.exit(signum)

    # store the original SIGTERM signal handler to restore it
    # once search systems is over
    default_signal_handler = signal.signal(signal.SIGTERM, stop)


    def search(gene, cpu):
        """
        Search gene in the database built from the input sequence file (execute \"hmmsearch\"), and produce a HMMReport

        :param gene: the gene to search
        :type gene: a :class:`macsypy.gene.CoreGene` object
        :param int cpu: the number of cpu to use (per worker)
        """
        _log.info(f"search gene {gene.name}")
        profile = gene.profile
        try:
            report = profile.execute(cpu=cpu)
        except Exception as err:
            _log.critical(err)
            stop(signal.SIGKILL, None)
        else:
            if report:
                report.extract()
                report.save_extract()
                return report


    def recover(gene, cfg):
        """
        Recover Hmmer output from a previous run, and produce a report

        :param gene: the gene to search
        :type gene: a :class:`macsypy.gene.CoreGene` object
        :param cfg: the configuration
        :type cfg: :class:`macsypy.config.Config` object
        :return: the list of all HMMReports (derived class depending on the input dataset type)
        :rtype: list of `macsypy.report.HMMReport` object
        """
        hmm_old_path = os.path.join(cfg.previous_run(), cfg.hmmer_dir(), gene.name + cfg.res_search_suffix())
        _log.info(f"recover hmm {hmm_old_path}")
        hmm_new_path = os.path.join(cfg.working_dir(), cfg.hmmer_dir(), gene.name + cfg.res_search_suffix())
        shutil.copy(hmm_old_path, hmm_new_path)
        gene.profile.hmm_raw_output = hmm_new_path
        db_type = cfg.db_type()
        if db_type == 'gembase':
            report = GembaseHMMReport(gene, hmm_new_path, cfg)
        elif db_type == 'ordered_replicon':
            report = OrderedHMMReport(gene, hmm_new_path, cfg)
        else:
            report = GeneralHMMReport(gene, hmm_new_path, cfg)
        if report:
            report.extract()
            report.save_extract()
            return report

    # there is only one instance of gene per name but the same instance can be
    # in all genes several times
    # hmmsearch and extract should be executed ONLY ONCE per run
    # so I uniquify the list of gene (use CoreGene and set)
    genes = {mg.core_gene for mg in genes}
    _log.debug("start searching genes")

    hmmer_dir = os.path.join(cfg.working_dir(), cfg.hmmer_dir())
    if not os.path.exists(hmmer_dir):
        # it works because mkdir is an atomic operation
        os.mkdir(hmmer_dir)

    # ensure that the indexes is build only once
    # So threads just have to read
    idx = Indexes(cfg)
    idx.build()
    previous_run = cfg.previous_run()
    with concurrent.futures.ThreadPoolExecutor(max_workers=worker_nb) as executor:
        future_search = []
        for gene in genes:
            if previous_run and os.path.exists(os.path.join(previous_run, cfg.hmmer_dir(),
                                                            gene.name + cfg.res_search_suffix())):
                future_search.append(executor.submit(recover, gene, cfg))
            else:
                future_search.append(executor.submit(search, gene, cpu_per_worker))
        for future in concurrent.futures.as_completed(future_search):
            report = future.result()
            if report:
                all_reports.append(report)
    _log.debug("end searching genes")

    signal.signal(signal.SIGTERM, default_signal_handler)
    return all_reports
