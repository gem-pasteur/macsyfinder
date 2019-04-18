# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur, Paris.                                   #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


from operator import attrgetter
from macsypy import model
from macsypy.hit import Hit
from macsypy.definition_parser import DefinitionParser
from macsypy.search_systems import Cluster
from macsypy.search_genes import search_genes
from macsypy.error import SystemDetectionError
from tests import MacsyTest
from tests.macsy_test_env import MacsyEnvManager


class Test(MacsyTest, MacsyEnvManager):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_len(self):
        self.load_env("env_001", log_out=False)

        system = model.Model(self.macsy_test_env.cfg, 'foo', 10)
        cluster = Cluster(system, self.macsy_test_env.model_bank)
        li = [None] * 12
        hit = Hit(*li)
        cluster.hits = [hit]
        self.assertEqual(len(cluster), 1)

        self.unload_env("env_001")

    def test_putative_system(self):
        self.load_env("env_001", log_out=False)

        model_name = 'set_1/T9SS'
        t9ss_model = model.Model(self.macsy_test_env.cfg, model_name, 10)
        cluster = Cluster(t9ss_model, self.macsy_test_env.model_bank)
        cluster._putative_system = model_name
        self.assertEqual(cluster.putative_system, model_name)

        self.unload_env("env_001")

    def test_compatible_systems(self):
        self.load_env("env_001", log_out=False)

        t9ss_model = model.Model(self.macsy_test_env.cfg, 'set_1/T9SS', 10)
        cluster = Cluster(t9ss_model, self.macsy_test_env.model_bank)
        compatible_system_name = 'set_1/T2SS'
        cluster._compatible_systems.append(compatible_system_name)
        self.assertListEqual(cluster.compatible_systems, [compatible_system_name])

        self.unload_env("env_001")

    def test_state(self):
        self.load_env("env_001", log_out=False)

        system = model.Model(self.macsy_test_env.cfg, 'foo', 4)
        cluster = Cluster(system, self.macsy_test_env.model_bank)
        state = cluster.state
        self.assertEqual(state, '')

        self.unload_env("env_001")

    def test_add(self):
        self.load_env("env_007", log_out=False)

        h1, h2, h3 = self.macsy_test_env.all_hits[:3]
        t9ss_model = self.macsy_test_env.models[0]
        cluster = Cluster([t9ss_model], self.macsy_test_env.model_bank)
        cluster.add(h1)
        self.assertEqual(cluster.begin, 505)
        self.assertEqual(cluster.end, 505)
        self.assertEqual(cluster.replicon_name, 'AESU001c01a')
        self.assertEqual(len(cluster.hits), 1)

        cluster.add(h2)
        self.assertEqual(cluster.begin, 505)
        self.assertEqual(cluster.end, 773)
        self.assertEqual(cluster.replicon_name, 'AESU001c01a')
        self.assertEqual(len(cluster.hits), 2)

        h3.replicon_name = 'bar'
        with self.assertRaises(SystemDetectionError) as context:
            cluster.add(h3)
        self.assertEqual(str(context.exception),
                         "Attempting to gather in a cluster hits from different replicons ! ")

        self.unload_env("env_007")

    def test_save(self):

        def get_hits(system):
            genes = system.mandatory_genes + system.accessory_genes + system.forbidden_genes

            ex_genes = []
            for g in genes:
                if g.exchangeable:
                    ex_genes += g.get_homologs() + g.get_analogs()
            all_genes = genes + ex_genes
            all_reports = search_genes(all_genes, self.macsy_test_env.cfg)
            all_hits = [hit for subl in [report.hits for report in all_reports] for hit in subl]
            all_hits = sorted(all_hits, key=attrgetter('score'), reverse=True)
            all_hits = sorted(all_hits, key=attrgetter('replicon_name', 'position'))

            return all_hits

        self.load_env("env_008", log_out=False)
        model_bank = self.macsy_test_env.model_bank
        gene_bank = self.macsy_test_env.gene_bank
        parser = DefinitionParser(self.macsy_test_env.cfg, model_bank, gene_bank,
                                  self.macsy_test_env.profile_factory)
        parser.parse(['set_1/T9SS'])
        model_1 = model_bank['set_1/T9SS']
        all_hits_1 = get_hits(model_1)

        # test case 1
        cluster = Cluster([model_1], model_bank)
        for h in all_hits_1:
            cluster.add(h)
        cluster.save()
        self.assertEqual(cluster.systems_count, {'set_1/T9SS': 17})
        self.assertEqual(cluster._putative_system, 'set_1/T9SS')
        self.assertEqual(cluster._compatible_systems, ['set_1/T9SS'])
        self.assertEqual(cluster._state, 'clear')

        # test case 2
        cluster = Cluster([model_1], self.macsy_test_env.model_bank)
        cluster.add(all_hits_1[0])
        cluster.hits[0].gene._loner = False
        cluster.save()
        self.assertEqual(cluster.systems_count, {'set_1/T9SS': 1})
        self.assertEqual(cluster._putative_system, 'set_1/T9SS')
        self.assertEqual(cluster._compatible_systems, [])
        self.assertEqual(cluster._state, 'ineligible')

        # test case 3
        fqn = 'set_1/T3SS'
        parser.parse([fqn])
        model_2 = model_bank[fqn]
        all_hits_2 = get_hits(model_2)
        cluster = Cluster([model_1, model_2], self.macsy_test_env.model_bank)
        for h in all_hits_1 + all_hits_2:
            cluster.add(h)
        cluster.save()
        self.assertEqual(cluster.systems_count, {'set_1/T9SS': 17, 'set_1/T3SS': 2})
        self.assertEqual(cluster._putative_system, 'set_1/T9SS')
        self.assertEqual(cluster._compatible_systems, [])
        self.assertEqual(cluster._state, 'ambiguous')

        self.unload_env("env_008")

        # test case 4
        self.load_env("env_006", log_out=False)
        model_bank = self.macsy_test_env.model_bank
        gene_bank = self.macsy_test_env.gene_bank
        parser = DefinitionParser(self.macsy_test_env.cfg, model_bank, gene_bank, self.macsy_test_env.profile_factory)
        fqn_1 = 'set_1/T2SS'
        fqn_2 = 'set_1/T4P'
        parser.parse([fqn_1, fqn_2])
        model_1 = model_bank[fqn_1]
        model_2 = model_bank[fqn_2]
        all_hits_1 = get_hits(model_1)
        all_hits_2 = get_hits(model_2)

        hits = {}
        for h in all_hits_1 + all_hits_2:
            if h.id in ['VICH001.B.00001.C001_00829', 'VICH001.B.00001.C001_00830', 'VICH001.B.00001.C001_00833']:
                if h.id not in hits:
                    hits[h.id] = h

        cluster = Cluster([model_1, model_2], self.macsy_test_env.model_bank)
        for h in hits.values():
            cluster.add(h)
        cluster.save()

        self.unload_env("env_006")

    def test_str(self):
        self.load_env("env_002", log_out=False)
        buffer_ = str(self.macsy_test_env.cluster)
        self.assertEqual(str(buffer_), self.output_control_str('001'))
        self.unload_env("env_002")
