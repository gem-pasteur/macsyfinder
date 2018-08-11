import os
import shutil
import tempfile
import logging
from time import strftime
from operator import attrgetter
from macsypy.database import Indexes, RepliconDB
from macsypy.config import Config
from macsypy.registries import ModelRegistry
from macsypy.gene import gene_bank
from macsypy.system import system_bank
from macsypy.search_systems import system_name_generator, build_clusters, analyze_clusters_replicon
from macsypy.search_genes import search_genes
from macsypy.system_parser import SystemParser
from tests import MacsyTest

class MacsyTestEnvSnippet():
    def build_hits(self):
        self.out_dir = "/tmp/macsyfinder-test_fill_with_cluster-" + strftime("%Y%m%d_%H-%M-%S")

        self.cfg = Config(hmmer_exe="hmmsearch",
                        out_dir=self.out_dir,
                        db_type="gembase",
                        previous_run="tests/data/data_set_1/complete_run_results",
                        e_value_res=1,
                        i_evalue_sel=0.5,
                        res_search_suffix=".search_hmm.out",
                        profile_suffix=".hmm",
                        res_extract_suffix="",
                        log_level=30,
                        models_dir="tests/data/data_set_1/models",
                        log_file=os.devnull)

        idx = Indexes(self.cfg)
        idx._build_my_indexes()

        parser = SystemParser(self.cfg, system_bank, gene_bank)
        parser.parse(['set_1/T9SS'])

        self.system = system_bank['set_1/T9SS']

        genes = self.system.mandatory_genes + self.system.accessory_genes + self.system.forbidden_genes

        ex_genes = []
        for g in genes:
            if g.exchangeable:
                h_s = g.get_homologs()
                ex_genes += h_s
                a_s = g.get_analogs()
                ex_genes += a_s
        all_genes = (genes + ex_genes)

        all_reports = search_genes(all_genes, self.cfg)

        all_hits = [hit for subl in [report.hits for report in all_reports] for hit in subl]

        all_hits = sorted(all_hits, key=attrgetter('score'), reverse=True)
        self.all_hits = sorted(all_hits, key=attrgetter('replicon_name', 'position'))

class MacsyTestEnv(MacsyTestEnvSnippet):
    """Standard test environments.
    
    env_001 => environment loaded from scratch 
               (data from "test_base.fa")
    env_002 => environment loaded using "previous_run" option
               (data from "tests/data/data_set_1/complete_run_results")
    """

    def load(self, env_id):
        l = logging.getLogger()
        l.manager.loggerDict.clear()

        # add only one handler to the macsypy logger
        from macsypy.report import _log
        macsy_log = _log.parent
        log_handler = logging.FileHandler(os.devnull)
        macsy_log.addHandler(log_handler)

        if env_id == "env_001":
            self.cfg = Config(hmmer_exe="hmmsearch",
                              sequence_db=MacsyTest.find_data("base", "test_base.fa"),
                              db_type="gembase",
                              e_value_res=1,
                              i_evalue_sel=0.5,
                              res_search_dir=tempfile.gettempdir(),
                              res_search_suffix=".search_hmm.out",
                              profile_suffix=".hmm",
                              res_extract_suffix="",
                              log_level=30,
                              models_dir=MacsyTest.find_data('models'),
                              log_file=os.devnull)

            shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
            self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))

            idx = Indexes(self.cfg)
            idx._build_my_indexes()

            models_registry = ModelRegistry(self.cfg)
            self.model_name = 'foo'
            self.models_location = models_registry[self.model_name]
        elif env_id == "env_002":
            self.build_hits()

            # debug
            # print [h.gene.name for h in self.all_hits]

            rep_db = RepliconDB(self.cfg)
            self.rep_info = rep_db['AESU001c01a']

            (clusters, multi_syst_genes) = build_clusters(self.all_hits, [self.system], self.rep_info)
            self.cluster = clusters.clusters[0]

            systems_occurences_list = analyze_clusters_replicon(clusters, [self.system], multi_syst_genes)

            self.system_occurence = systems_occurences_list[0]

            # debug
            #print system_occurence.valid_hits
        elif env_id == "env_003":
            self.build_hits()

            # debug
            # print [h.gene.name for h in self.all_hits]

            rep_db = RepliconDB(self.cfg)
            self.rep_info = rep_db['AESU001c01a']
        else:
            raise Exception('Test environment not found ({})'.format(env_id))

    def unload(self, env_id):

        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

        # reset global vars
        system_bank._system_bank = {}
        gene_bank._genes_bank = {}
        system_name_generator.name_bank = {}

        # environment specific cleanup
        if env_id == "env_001":
            pass
        elif env_id == "env_002":
            try:
                shutil.rmtree(self.out_dir)
            except:
                pass
        elif env_id == "env_003":
            try:
                shutil.rmtree(self.out_dir)
            except:
                pass
        else:
            raise Exception('Test environment not found ({})'.format(env_id))
