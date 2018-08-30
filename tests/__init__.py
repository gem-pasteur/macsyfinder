import os
import sys
import shutil
import unittest
import platform
from StringIO import StringIO
from contextlib import contextmanager
import hashlib
from functools import partial
import tempfile
import uuid
import logging
from operator import attrgetter

from macsypy.database import Indexes, RepliconDB
from macsypy.config import Config
from macsypy.registries import ModelRegistry
from macsypy.gene import gene_bank
from macsypy.system import system_bank
from macsypy.search_systems import system_name_generator, build_clusters, analyze_clusters_replicon
from macsypy.search_genes import search_genes
from macsypy.system_parser import SystemParser


class MacsyTest(unittest.TestCase):

    _tests_dir = os.path.normpath(os.path.dirname(__file__))
    _data_dir = os.path.join(_tests_dir, "data")
    _output_control_dir = os.path.join(_data_dir, "outputs_control")

    def setsid(self):
        platform = sys.platform
        if platform.startswith('linux'):
            setsid = 'setsid'
        elif platform.startswith('darwin'):
            setsid = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'utils', 'setsid'))
        else:
            setsid = ''
        return setsid

    @classmethod
    def find_data(cls, *args):
        data_path = os.path.join(cls._data_dir, *args)
        if os.path.exists(data_path):
            return data_path
        else:
            raise IOError("data '{}' does not exists".format(data_path))

    @contextmanager
    def catch_io(self, out=False, err=False):
        """
        Catch stderr and stdout of the code running within this block.
        """
        old_out = sys.stdout
        new_out = old_out
        old_err = sys.stderr
        new_err = old_err
        if out:
            new_out = StringIO()
        if err:
            new_err = StringIO()
        try:
            sys.stdout, sys.stderr = new_out, new_err
            yield sys.stdout, sys.stderr
        finally:
            sys.stdout, sys.stderr = old_out, old_err

    @staticmethod
    def fake_exit(*args, **kwargs):
        returncode = args[0]
        raise TypeError(returncode)

    @staticmethod
    def mute_call(call_ori):
        """
        hmmsearch or prodigal write lot of things on stderr or stdout
        which noise the unit test output
        So I replace the `call` function in module integron_finder
        by a wrapper which call the original function but add redirect stderr and stdout
        in dev_null
        :return: wrapper around call function
        :rtype: function
        """
        def wrapper(*args, **kwargs):
            with open(os.devnull, 'w') as f:
                kwargs['stderr'] = f
                kwargs['stdout'] = f
                res = call_ori(*args, **kwargs)
            return res
        return wrapper

    def assertFileEqual(self, f1, f2, msg=None):
        self.maxDiff = None
        # the StringIO does not support context in python2.7
        # so we can use the following statement only in python3
        # with open(f1) if isinstance(f1, str) else f1 as fh1, open(f2) if isinstance(f2, str) else f2 as fh2:
        with open(f1) as fh1, open(f2) as fh2:
            self.assertMultiLineEqual(fh1.read(), fh2.read(), msg=msg)

    def assertSeqRecordEqual(self, s1, s2):
        for attr in ('id', 'name', 'seq'):
            s1_attr = getattr(s1, attr)
            s2_attr = getattr(s2, attr)
            self.assertEqual(s1_attr, s2_attr, msg="{} are different: {} != {}".format(attr, s1_attr, s2_attr))

        # there is a bug in some biopython version
        self.assertEqual(s1.description.rstrip('.'), s2.description.rstrip('.'))
        for s1_feat, s2_feat in zip(s1.features, s2.features):
            # location cannot be directly compared
            self.assertEqual(str(s1_feat.location), str(s2_feat.location))

            for attr in ('qualifiers', 'strand', 'type'):
                f1_attr = getattr(s1_feat, attr)
                f2_attr = getattr(s2_feat, attr)
                self.assertEqual(f1_attr, f2_attr, msg="{} are different: {} != {}".format(attr, f1_attr, f2_attr))

    def assertHmmEqual(self, hmm1, hmm2):
        with open(hmm1) as hmm1_file, open(hmm2) as hmm2_file:
            for hmm1_line, hmm2_line in zip(hmm1_file, hmm2_file):
                if hmm1_line.startswith('#') and hmm2_line.startswith('#'):
                    continue
                hmm1_fields = hmm1_line.split('#')[:-1]
                hmm2_fields = hmm2_line.split('#')[:-1]
                self.assertListEqual(hmm1_fields, hmm2_fields)

    @staticmethod
    def get_uniq_tmp_dir_name():
        return os.path.join(tempfile.gettempdir(), "macsyfinder-{}".format(uuid.uuid4()))

    @staticmethod
    def rmtree(path):
        """
        Remove directory tree.

        :param path: the path to remove
        :type path: str
        """
        try:
            shutil.rmtree(path)
        except:
            pass

    @staticmethod
    def md5sum(file_=None, str_=None):
        """Compute md5 checksum.

        :param file_: the name of the file to compute the checksum for
        :type file_: str
        :param str_: the string to compute the checksum for
        :type str_: str
        """
        assert not (file_ and str_)

        d = hashlib.md5()

        if file_:
            with open(file_, mode='rb') as f:
                for buf in iter(partial(f.read, 128), b''):
                    d.update(buf)
        elif str_:
            assert isinstance(str_, str)
            d.update(str_)
        else:
            assert False

        return d.hexdigest()


class LoggerWrapper(object):

    def __init__(self, logger):
        self.logger = logger

    def __getattr__(self, item):
        return getattr(self.logger, item)

    def get_value(self):
        return self.logger.handlers[0].stream.getvalue()


class MacsyTestEnvSnippet():

    def build_config(self, previous_run="tests/data/data_set_3/results", models_dir="tests/data/data_set_3/models"):
        self.out_dir = MacsyTest.get_uniq_tmp_dir_name()

        self.config = Config(hmmer_exe="hmmsearch",
                             out_dir=self.out_dir,
                             db_type="gembase",
                             previous_run=previous_run,
                             e_value_res=1,
                             i_evalue_sel=0.5,
                             res_search_suffix=".search_hmm.out",
                             profile_suffix=".hmm",
                             res_extract_suffix="",
                             log_level=30,
                             models_dir=models_dir,
                             log_file=os.devnull)

        idx = Indexes(self.config)
        idx._build_my_indexes()

    def build_hits(self, previous_run="tests/data/data_set_1/complete_run_results", models_dir="tests/data/data_set_1/models"):
        self.out_dir = MacsyTest.get_uniq_tmp_dir_name()

        self.cfg = Config(hmmer_exe="hmmsearch",
                          out_dir=self.out_dir,
                          db_type="gembase",
                          previous_run=previous_run,
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          res_search_suffix=".search_hmm.out",
                          profile_suffix=".hmm",
                          res_extract_suffix="",
                          log_level=30,
                          models_dir=models_dir,
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
            # print system_occurence.valid_hits
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
            self.build_hits(previous_run="tests/data/data_set_2/results", models_dir="tests/data/data_set_2/models")
        elif env_id == "env_006":
            self.build_config(previous_run="tests/data/data_set_3/results", models_dir="tests/data/data_set_3/models")
        elif env_id == "env_007":
            self.build_hits()
        elif env_id == "env_008":
            self.build_config(previous_run="tests/data/data_set_1/complete_run_results", models_dir="tests/data/data_set_1/models")
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
        else:
            raise Exception('Test environment not found ({})'.format(env_id))


def which(name, flags=os.X_OK):
    """
    Search PATH for executable files with the given name.

    :param name: the name of the executable to search
    :type name: str
    :param flags: os mod the name must have, default is executable (os.X_OK).
    :type flags: os file mode R_OK|R_OK|W_OK|X_OK
    :return: the path of the executable
    :rtype: string or None
    """
    result = None
    path = os.environ.get('PATH', None)
    if path is None:
        return result
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if platform.system() == 'Windows':
            p += '.exe'
        if os.access(p, flags):
            result = p
            break
    return result
