import os
import shutil
import tempfile
import logging
from macsypy.database import Indexes
from macsypy.config import Config
from macsypy.registries import ModelRegistry
from macsypy.gene import gene_bank
from macsypy.system import system_bank
from macsypy.search_systems import system_name_generator
from tests import MacsyTest

class MacsyTestEnv():
    """Standard test environments."""

    def load(self, env_id):
        if env_id == "env_001":
            l = logging.getLogger()
            l.manager.loggerDict.clear()

            # add only one handler to the macsypy logger
            from macsypy.report import _log
            macsy_log = _log.parent
            log_file = os.devnull
            log_handler = logging.FileHandler(log_file)
            macsy_log.addHandler(log_handler)

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
                              log_file=log_file)

            shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
            self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))

            idx = Indexes(self.cfg)
            idx._build_my_indexes()

            models_registry = ModelRegistry(self.cfg)
            self.model_name = 'foo'
            self.models_location = models_registry[self.model_name]
        else:
            raise Exception('Test environment not found ({})'.format(env_id))

    def unload(self, env_id):
        if env_id == "env_001":
            # close loggers filehandles, so they don't block file deletion
            # in shutil.rmtree calls in Windows
            logging.shutdown()
            l = logging.getLogger()
            l.manager.loggerDict.clear()
            try:
                shutil.rmtree(self.cfg.working_dir)
            except:
                pass
        else:
            raise Exception('Test environment not found ({})'.format(env_id))

        system_bank._system_bank = {}
        gene_bank._genes_bank = {}
        system_name_generator.name_bank = {}
