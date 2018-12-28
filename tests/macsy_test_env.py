import os
import shutil
import tempfile
import logging
from operator import attrgetter
import argparse

from macsypy.database import Indexes, RepliconDB
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelRegistry
from macsypy.gene import gene_bank
from macsypy.system import system_bank
from macsypy.search_systems import system_name_generator, build_clusters, analyze_clusters_replicon
from macsypy.search_genes import search_genes
from macsypy.system_parser import SystemParser

from tests import MacsyTest


class MacsyTestEnvSnippet(object):

    def __init__(self):
        self.out_dir = None
        self.cfg = None
        self.system = None
        self.all_hits = []

    def build_config(self, previous_run="tests/data/data_set_3/results",
                           models_dir="tests/data/data_set_3/models",
                           sequence_db="test_base.fa"):
        self.out_dir = MacsyTest.get_uniq_tmp_dir_name()
        seq_ori = MacsyTest.find_data("base", sequence_db)

        defaults = MacsyDefaults()
        args = argparse.Namespace()
        args.out_dir = self.out_dir
        args.db_type = 'gembase'
        args.previous_run = previous_run
        args.log_level = 30
        args.models_dir = models_dir
        args.log_file = os.devnull
        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        os.mkdir(args.out_dir)
        args.sequence_db = os.path.join(args.out_dir, os.path.basename(seq_ori))
        shutil.copy(seq_ori, args.out_dir)
        self.cfg = Config(defaults, args)
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~ previous_run", self.cfg.previous_run())
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~ sequence_db", self.cfg.sequence_db())
        idx = Indexes(self.cfg)
        idx.build()


    def build_hits(self, previous_run="tests/data/data_set_1/complete_run_results",
                   models_dir="tests/data/data_set_1/models",
                   model_fqn="set_1/T9SS"):
        self.build_config(previous_run=previous_run, models_dir=models_dir)
        parser = SystemParser(self.cfg, system_bank, gene_bank)
        parser.parse([model_fqn])
        self.system = system_bank[model_fqn]
        genes = self.system.mandatory_genes + self.system.accessory_genes + self.system.forbidden_genes

        ex_genes = []
        for g in genes:
            if g.exchangeable:
                h_s = g.get_homologs()
                ex_genes += h_s
                a_s = g.get_analogs()
                ex_genes += a_s
        all_genes = (genes + ex_genes)
        print("\n########## all_genes", [g.name for g in all_genes] )
        all_reports = search_genes(all_genes, self.cfg)
        print("\n################# all_reports", all_reports)
        all_hits = [hit for subl in [report.hits for report in all_reports] for hit in subl]
        all_hits = sorted(all_hits, key=attrgetter('score'), reverse=True)
        self.all_hits = sorted(all_hits, key=attrgetter('replicon_name', 'position'))


class MacsyTestEnv(MacsyTestEnvSnippet):
    """Standard test environments.

    env_001 => environment loaded from scratch
               (data from "test_base.fa").
    env_002 => environment loaded using "previous_run" option
               (data from "tests/data/data_set_1/complete_run_results").
               Create SystemOccurence.
    env_003 => environment loaded using "previous_run" option
               (data from "tests/data/data_set_1/complete_run_results").
               Stops before calling build_clusters() method.
    env_004 => environment loaded using "previous_run" option
               (data from "tests/data/data_set_1/complete_run_results").
               Create ModelRegistry.
    env_005 => environment loaded using "previous_run" option
               (data from "tests/data/data_set_2/results").
               Stops before calling build_clusters() method.
               Do not contain RepliconDB instance.
    env_006 => environment loaded using "previous_run" option
               (use data from "tests/data/data_set_3").
               Stops after creating indexes.
    env_007 => environment loaded using "previous_run" option.
               (use data from "tests/data/data_set_1").
               Stops after creating hits.
    env_008 => environment loaded using "previous_run" option
               (use data from "tests/data/data_set_1").
               Stops after creating indexes.
    env_009 => environment loaded using "previous_run" option
               (use data from "tests/data/data_set_3").
               Stops after creating ModelRegistry.
    env_010 => environment loaded from scratch
               (data from "test_base_with_errors.fa").
               Index not created.
    """


    def __init__(self):
        super(MacsyTestEnv, self).__init__()
        self.model_name = None
        self.models_location = None
        self.rep_info = None
        self.cluster = None
        self.system_occurence = None
        self.defaults = MacsyDefaults()

    def load(self, env_id):
        l = logging.getLogger()
        logging._handlers.clear()                                                                                                                                                                                             
        logging.shutdown(logging._handlerList[:])                                                                                                                                                                             
        del logging._handlerList[:]                                                                                                                                                                                           
        l.manager.loggerDict.clear()

        # add only one handler to the macsypy logger
        from macsypy.report import _log
        macsy_log = _log.parent
        log_handler = logging.FileHandler(os.devnull)
        macsy_log.addHandler(log_handler)

        if env_id == "env_001":
            self.build_config(previous_run=None,
                              models_dir=MacsyTest.find_data('models')
                              )
            idx = Indexes(self.cfg)
            idx.build()

            models_registry = ModelRegistry(self.cfg)
            self.model_name = 'foo'
            self.models_location = models_registry[self.model_name]

        elif env_id == "env_002":
            self.build_hits()
            rep_db = RepliconDB(self.cfg)
            self.rep_info = rep_db['AESU001c01a']
            clusters, multi_syst_genes = build_clusters(self.all_hits, [self.system], self.rep_info)
            self.cluster = clusters.clusters[0]
            systems_occurences_list = analyze_clusters_replicon(clusters, [self.system], multi_syst_genes)
            self.system_occurence = systems_occurences_list[0]

        elif env_id == "env_003":
            self.build_hits()
            rep_db = RepliconDB(self.cfg)
            self.rep_info = rep_db['AESU001c01a']
        elif env_id == "env_004":
            self.build_hits()
            models_registry = ModelRegistry(self.cfg)
            self.model_name = 'set_1'
            self.models_location = models_registry[self.model_name]
        elif env_id == "env_005":
            self.build_hits(previous_run="tests/data/data_set_2/results",
                            models_dir="tests/data/data_set_2/models")
        elif env_id == "env_006":
            self.build_config(previous_run="tests/data/data_set_3/results",
                              models_dir="tests/data/data_set_3/models")
        elif env_id == "env_007":
            self.build_hits()
        elif env_id == "env_008":
            self.build_config(previous_run="tests/data/data_set_1/complete_run_results",
                              models_dir="tests/data/data_set_1/models")
        elif env_id == "env_009":
            self.build_hits(previous_run="tests/data/data_set_3/results",
                            models_dir="tests/data/data_set_3/models",
                            model_fqn="set_1/T4P")

            models_registry = ModelRegistry(self.cfg)
            self.model_name = 'set_1'
            self.models_location = models_registry[self.model_name]
        elif env_id == "env_010":
            args = argparse.Namespace()
            args.sequence_db = MacsyTest.find_data("base",  "test_base_with_errors.fa")
            args.db_type = 'gembase'
            args.res_search_dir = MacsyTest.get_uniq_tmp_dir_name()
            args.log_level = 30
            args.log_file = os.devnull
            args.models_dir = MacsyTest.find_data('models')
            self.cfg = Config(self.defaults, args)
        else:
            raise Exception('Test environment not found ({})'.format(env_id))

    def unload(self, env_id):

        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        MacsyTest.close_loggers_filehandles()
        try:
            shutil.rmtree(self.cfg.working_dir)
        except Exception:
            pass

        # reset global vars
        system_bank._system_bank = {}
        gene_bank._genes_bank = {}
        system_name_generator.name_bank = {}

        # environment specific cleanup
        if env_id == "env_001":
            pass
        elif env_id == "env_002":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_003":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_004":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_005":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_006":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_007":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_008":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_009":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_010":
            pass
        else:
            raise Exception('Test environment not found ({})'.format(env_id))


class MacsyEnvManager(object):

    def __init__(self):
        self.macsy_test_env = None

    def load_env(self, env_id):
        self.macsy_test_env = MacsyTestEnv()
        self.macsy_test_env.load(env_id)

    def unload_env(self, env_id):
        self.macsy_test_env.unload(env_id)
        self.macsy_test_env = None
