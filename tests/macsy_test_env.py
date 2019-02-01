import os
import shutil
import logging
from operator import attrgetter
import argparse

from macsypy.database import Indexes, RepliconDB
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelRegistry
from macsypy import gene
from macsypy import model
from macsypy.search_systems import build_clusters, analyze_clusters_replicon
from macsypy import search_systems
from macsypy.search_genes import search_genes
from macsypy.definition_parser import DefinitionParser
from macsypy.utils import get_models_name_to_detect
from tests import MacsyTest
import macsypy


class MacsyTestEnvSnippet(object):

    def __init__(self):
        self.out_dir = None
        self.cfg = None
        self.models = []
        self.all_hits = []
        self.gene_bank = gene.GeneBank()
        gene.gene_bank = self.gene_bank

        self.model_bank = model.ModelBank()
        model.model_bank = self.model_bank
        search_systems.model_bank = self.model_bank

    def build_config(self, **config_opts):
        assert self.out_dir is not None

        defaults = MacsyDefaults()
        seq_ori = MacsyTest.find_data("base", config_opts.get('sequence_db', 'test_base.fa'))
        if 'sequence_db' in config_opts:
            del config_opts['sequence_db']

        args = argparse.Namespace()
        args.out_dir = self.out_dir
        args.db_type = 'gembase'
        args.log_level = 20
        args.models_dir = MacsyTest.find_data("data_set_3", "models")
        args.log_file = os.devnull
        for opt, v in config_opts.items():
            setattr(args, opt, v)
        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        os.mkdir(args.out_dir)
        if 'previous_run' not in config_opts:
            args.sequence_db = os.path.join(args.out_dir, os.path.basename(seq_ori))
            shutil.copy(seq_ori, args.out_dir)
        self.cfg = Config(defaults, args)
        idx = Indexes(self.cfg)
        idx.build()
        self.registry = ModelRegistry(self.cfg)

    def build_hits(self, models_2_parse=None, **config_opts):
        default_opts = {'previous_run': "tests/data/data_set_1/complete_run_results",
                        'models_dir': "tests/data/data_set_1/models"
                        }
        default_opts.update(config_opts)
        self.build_config(**default_opts)
        parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank)
        if models_2_parse is None:
            models_2_parse = get_models_name_to_detect(self.cfg.models(), self.registry)
        parser.parse(models_2_parse)
        self.models = [self.model_bank[model_fqn] for model_fqn in models_2_parse]
        genes = []
        for model in self.models:
            genes += model.mandatory_genes + model.accessory_genes + model.forbidden_genes
        genes = list(set(genes))
        ex_genes = []
        for g in genes:
            if g.exchangeable:
                ex_genes += g.get_homologs() + g.get_analogs()
        all_genes = genes + ex_genes
        all_reports = search_genes(all_genes, self.cfg)
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
    env_011 => environment loaded using "previous_run" option
               (data from "tests/data/data_set_4/results").
               Stops before calling build_clusters() method.
               Do not contain RepliconDB instance.
    env_012 => environment loaded using "previous_run" option
               (data from "tests/data/data_set_5/results").
               Stops before calling build_clusters() method.
               Do not contain RepliconDB instance.
    env_013 => environment loaded using "previous_run" option
               (data from "tests/data/data_set_6/results").
               Stops before calling build_clusters() method.
               Do not contain RepliconDB instance.
    """


    def __init__(self):
        super(MacsyTestEnv, self).__init__()
        self.model_name = None
        self.models_location = None
        self.rep_info = None
        self.cluster = None
        self.system_occurence = None
        self.defaults = MacsyDefaults()
        self.handlers = []


    def load(self, env_id, log_out=True, log_level=logging.INFO, **cfg_args):
        self.out_dir = MacsyTest.get_tmp_dir_name()

        MacsyTest.rmtree(self.out_dir)

        self.handlers = macsypy.init_logger(out=log_out)
        macsypy.logger_set_level(level=log_level)

        if env_id == "env_001":
            cfg_args.update({'models_dir': MacsyTest.find_data('models')})
            self.build_config(**cfg_args)
            idx = Indexes(self.cfg)
            idx.build()

            models_registry = ModelRegistry(self.cfg)
            self.model_name = 'foo'
            self.models_location = models_registry[self.model_name]

        elif env_id == "env_002":
            self.build_hits(previous_run='tests/data/data_set_1/complete_run_results',
                            **cfg_args)
            rep_db = RepliconDB(self.cfg)
            self.rep_info = rep_db['AESU001c01a']
            clusters, multi_syst_genes = build_clusters(self.all_hits, self.models, self.rep_info)
            self.cluster = clusters.clusters[0]
            systems_occurences_list = analyze_clusters_replicon(clusters, self.models, multi_syst_genes)
            self.system_occurence = systems_occurences_list[0]

        elif env_id == "env_003":
            self.build_hits(previous_run='tests/data/data_set_1/complete_run_results',
                            models_2_parse=["set_1/T3SS", "set_1/T4SS_typeI", "set_1/T9SS"],
                            **cfg_args)
            rep_db = RepliconDB(self.cfg)
            self.rep_info = rep_db['AESU001c01a']
        elif env_id == "env_004":
            self.build_hits(previous_run='tests/data/data_set_1/complete_run_results',
                            **cfg_args)
            models_registry = ModelRegistry(self.cfg)
            self.model_name = 'set_1'
            self.models_location = models_registry[self.model_name]
        elif env_id == "env_005":
            self.build_hits(previous_run="tests/data/data_set_2/results",
                            models_dir="tests/data/data_set_2/models",
                            models_2_parse=["set_1/T9SS"],
                            **cfg_args)
        elif env_id == "env_006":
            self.build_config(previous_run="tests/data/data_set_3/results",
                              models_dir="tests/data/data_set_3/models",
                              model_fqn=["set_1/T9SS"],
                              **cfg_args)
        elif env_id == "env_007":
            self.build_hits(previous_run="tests/data/data_set_1/complete_run_results",
                            **cfg_args)
        elif env_id == "env_008":
            self.build_config(previous_run="tests/data/data_set_1/complete_run_results",
                              models_dir="tests/data/data_set_1/models",
                              **cfg_args)
        elif env_id == "env_009":
            self.build_hits(previous_run="tests/data/data_set_3/results",
                            models_dir="tests/data/data_set_3/models",
                            models_2_parse=["set_1/T4P"],
                            i_evalue_sel=0.5,
                            **cfg_args)

            models_registry = ModelRegistry(self.cfg)
            self.model_name = 'set_1'
            self.models_location = models_registry[self.model_name]
        elif env_id == "env_010":
            args = argparse.Namespace()
            args.sequence_db = MacsyTest.find_data("base",  "test_base_with_errors.fa")
            args.db_type = 'gembase'
            args.res_search_dir = self.out_dir
            args.log_level = 30
            args.log_file = os.devnull
            args.models_dir = MacsyTest.find_data('models')
            self.cfg = Config(self.defaults, args)
        elif env_id == "env_011":
            self.build_hits(previous_run="tests/data/data_set_4/results",
                            models_dir="tests/data/data_set_4/models",
                            models_2_parse=["set_1/T9SS"],
                            **cfg_args)
        elif env_id == "env_012":
            self.build_hits(previous_run="tests/data/data_set_5/results",
                            models_dir="tests/data/data_set_5/models",
                            models_2_parse=["set_1/T9SS"],
                            **cfg_args)
        elif env_id == "env_013":
            self.build_hits(previous_run="tests/data/data_set_6/results",
                            models_dir="tests/data/data_set_6/models",
                            models_2_parse=["set_1/T9SS"],
                            **cfg_args)
        else:
            raise Exception('Test environment not found ({})'.format(env_id))

    def unload(self, env_id):
        # multiple call to init_logger will add handlers
        # so we need to clean handlers added
        logger = logging.getLogger('macsypy')
        for h in self.handlers:
            logger.removeHandler(h)

        MacsyTest.rmtree(self.cfg.working_dir())

        # reset global vars
        model.model_bank = model.ModelBank()
        gene.gene_bank = gene.GeneBank()
        search_systems.system_name_generator = search_systems.SystemNameGenerator()

        # environment specific cleanup
        if env_id == "env_001":
            MacsyTest.rmtree(self.out_dir)
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
        elif env_id == "env_011":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_012":
            MacsyTest.rmtree(self.out_dir)
        elif env_id == "env_013":
            MacsyTest.rmtree(self.out_dir)
        else:
            raise Exception('Test environment not found ({})'.format(env_id))


class MacsyEnvManager(object):

    def __init__(self):
        self.macsy_test_env = None

    def load_env(self, env_id, **cfg_args):
        self.macsy_test_env = MacsyTestEnv()
        self.macsy_test_env.load(env_id, **cfg_args)

    def unload_env(self, env_id):
        self.macsy_test_env.unload(env_id)
        self.macsy_test_env = None
