# -*- coding: utf-8 -*-

###############################
# Created on Jan 9, 2013
#
# @author: bneron
# @contact: bneron@pasteur.fr
# @organization: Institut Pasteur
# @license: GPLv3
################################

import threading


def search_genes(genes, cfg):
    """
    For each each gene of the list use the profile to perform an HMM and parse the output
    to generate a HMMReport and save it in a file. These tasks are performed in parrallel using threads.
    The number of workers can be limited by worker_nb directive in the config object.
    
    :param genes: the genes to search in the genome
    :type genes: list of :class:`txsscanlib.gene.Gene` objects
    :param cfg: the configuration 
    :type cfg: :class:`txsscanlib.config.Config` object
    """
    worker_nb = cfg.worker_nb
    if not worker_nb:
        worker_nb = len(genes)
    print "worker_nb =", worker_nb
    sema = threading.BoundedSemaphore(value = worker_nb)

    def worker(gene, sema):
        with sema:
            profile = gene.profile
            report = profile.execute()
            report.extract()
            report.save_extract()

    for g in genes:
        t = threading.Thread(target = worker, args = (g, sema))
        t.start()
    main_thread = threading.currentThread()    
    for t in threading.enumerate():
        if t is main_thread:
            continue
        t.join()
            