##########################################################################
#  MacSyFinder - Detection of macromolecular systems in protein dataset  #
#                using systems modelling and similarity search.          #
#  Authors: Sophie Abby, Bertrand Neron                                  #
#  Copyright (c) 2014-2022  Institut Pasteur (Paris) and CNRS.           #
#  See the COPYRIGHT file for details                                    #
#                                                                        #
#  This file is part of MacSyFinder package.                             #
#                                                                        #
#  MacSyFinder is free software: you can redistribute it and/or modify   #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  MacSyFinder is distributed in the hope that it will be useful,        #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details .                         #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with MacSyFinder (COPYING).                                     #
#  If not, see <https://www.gnu.org/licenses/>.                          #
##########################################################################
import argparse
import os.path
import io
import sys
import tempfile
import shutil
from unittest.mock import patch

from tests import MacsyTest

from macsypy import __version__ as msf_vers
from macsypy.config import MacsyDefaults, Config
from macsypy.error import MacsypyError
import macsypy.scripts.macsyconfig as msf_cfg


class TestConfigParserWithComment(MacsyTest):

    def test_add_comment(self):
        cp = msf_cfg.ConfigParserWithComments()
        cp.add_section('section_1')
        cp.add_comment('section_1', 'opt_1', 'comment_1')
        self.assertTrue(cp.has_option('section_1', 'opt_1_1_comment'))
        cp.add_comment('section_1', 'opt_1', 'comment_1')
        self.assertTrue(cp.has_option('section_1', 'opt_1_2_comment'))


    def test_write(self):

        cp = msf_cfg.ConfigParserWithComments()
        cp.add_section('section_1')
        cp.add_comment('section_1', 'opt_1', 'comment_1\ncomment_1 line 2', add_space_before=False, add_space_after=False)
        cp.add_comment('section_1', 'opt_1', 'comment_1_2', add_space_before=False, add_space_after=False)
        cp.set('section_1', 'opt_1', 'value_1')
        cp.add_comment('section_1', 'opt_2', 'comment_2', add_space_before=True, add_space_after=False)
        cp.set('section_1', 'opt_2', 'value_2')
        cp.add_comment('section_1', 'opt_3', 'comment_3', add_space_before=True, add_space_after=True)
        cp.set('section_1', 'opt_3', 'value_3')

        expected = """[section_1]
# comment_1
# comment_1 line 2
# comment_1_2
opt_1 = value_1

# comment_2
opt_2 = value_2

# comment_3

opt_3 = value_3

"""
        f = io.StringIO()
        cp.write(f)
        self.assertEqual(f.getvalue(),
                         expected)


class TestMacsyconfig(MacsyTest):

    def test_check_exe(self):
        exe = msf_cfg.check_exe('/bin/sh', None, None)
        self.assertEqual(exe, '/bin/sh')

        exe = msf_cfg.check_exe('', '/bin/sh', None)
        self.assertEqual(exe, '/bin/sh')

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_exe('nimportnaoik', None, None)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: 'nimportnaoik' NO executable found")


    def test_check_positive_int(self):
        val = msf_cfg.check_positive_int('1', None, None)
        self.assertEqual(val, 1)

        val = msf_cfg.check_positive_int('', 2, None)
        self.assertEqual(val, 2)

        val = msf_cfg.check_positive_int('1', None, None)
        self.assertEqual(val, 1)

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_positive_int('', None, None)
        self.assertEqual(str(ctx.exception),
                         'Please enter some value')

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_positive_int('-2', None, None)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: '-2' is not >=0")

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_positive_int('2.3', None, None)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: invalid literal for int() with base 10: '2.3'")

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_positive_int('nimportnaoik', None, None)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: invalid literal for int() with base 10: 'nimportnaoik'")


    def test_check_float(self):
        # send a good value
        val = msf_cfg.check_float('1.0', None, None)
        self.assertEqual(val, 1.0)

        # no value with default
        val = msf_cfg.check_float('', 1.0, None)
        self.assertEqual(val, 1.0)

        # check if casting goes well
        val = msf_cfg.check_float('1', None, None)
        self.assertEqual(val, 1.0)

        val = msf_cfg.check_float('', 2.1, None)
        self.assertEqual(val, 2.1)

        # # series of values
        val = msf_cfg.check_float('2.1, 3.2', '', None, sequence=True)
        self.assertEqual(val, [2.1, 3.2])

        # series of default values
        val = msf_cfg.check_float('', [2.1, 3.2], None, sequence=True)
        self.assertEqual(val, [2.1, 3.2])

        # no value, no default
        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_float('', None, None)
        self.assertEqual(str(ctx.exception),
                         'Please enter some value')
        # bad value
        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_float('nimportnaoik', None, None)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: could not convert string to float: 'nimportnaoik'")


    def test_check_str(self):
        val = msf_cfg.check_str("str", None, None)
        self.assertEqual(val, "str")

        # not declared as accepting sequence of string
        val = msf_cfg.check_str("str1, str2", None, None)
        self.assertEqual(val, "str1, str2")

        # not declared as accepting sequence of string
        val = msf_cfg.check_str("str1, str2", None, None, sequence=True)
        self.assertEqual(val, ["str1", "str2"])


    def test_check_bool(self):
        val = msf_cfg.check_bool('0', None, None)
        self.assertFalse(val)
        val = msf_cfg.check_bool('False', None, None)
        self.assertFalse(val)
        val = msf_cfg.check_bool('No', None, None)
        self.assertFalse(val)
        val = msf_cfg.check_bool('', 'No', None)
        self.assertFalse(val)

        val = msf_cfg.check_bool('1', None, None)
        self.assertTrue(val)
        val = msf_cfg.check_bool('True', None, None)
        self.assertTrue(val)
        val = msf_cfg.check_bool('Yes', None, None)
        self.assertTrue(val)
        val = msf_cfg.check_bool('', 'Yes', None)
        self.assertTrue(val)

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_bool('nimportnaoik', None, None)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: Authorized values ['True'/False/0/1]")


    def test_check_dir(self):
        test_dir = os.path.dirname(__file__)
        val = msf_cfg.check_dir(test_dir, None, None)
        self.assertEqual(val, test_dir)

        test_dir = os.path.dirname(__file__)
        val = msf_cfg.check_dir('', test_dir, None)
        self.assertEqual(val, test_dir)

        test_dir = __file__
        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_dir(test_dir, None, None)
        self.assertEqual(str(ctx.exception),
                         f"Invalid value: '{test_dir}' is not a directory.")

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_dir('nimportnaoik', None, None)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: 'nimportnaoik' no such file or directory.")

        test_dir = os.path.dirname(__file__)
        parent = os.path.join(os.path.normpath(os.path.join(test_dir, '..')))
        val = msf_cfg.check_dir(f"{test_dir}, {parent}", '', None, sequence=True)
        self.assertEqual(val, [test_dir, parent])


    def test_check_file(self):
        test_file = __file__
        sibling = os.path.join(os.path.dirname(test_file), 'test_hit.py')
        val = msf_cfg.check_file(test_file, None, None)
        self.assertEqual(val, test_file)

        val = msf_cfg.check_file('', test_file, None)
        self.assertEqual(val, test_file)

        val = msf_cfg.check_file('none', None, None)
        self.assertIsNone(val)

        val = msf_cfg.check_file(f"{test_file}, {sibling}", None, None, sequence=True)
        self.assertEqual(val, [test_file, sibling])

        test_file = os.path.dirname(__file__)
        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_file(test_file, None, None)
        self.assertEqual(str(ctx.exception),
                         f"Invalid value: '{test_file}' is not a file.")

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_file('nimportnaoik', None, None)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: 'nimportnaoik' no such file or directory.")


    def test_check_choice(self):
        choice = ['in', 'Yes']
        val = msf_cfg.check_choice('in', None, choice)
        self.assertEqual(val, 'in')

        val = msf_cfg.check_choice('', 'in', choice)
        self.assertEqual(val, 'in')

        with self.assertRaises(MacsypyError) as ctx:
            msf_cfg.check_choice('Out', None, choice)
        self.assertEqual(str(ctx.exception),
                         "Invalid value: Authorized values are ['in', 'Yes'].")

        val = msf_cfg.check_choice('none', None, ['in', 'none'])
        self.assertIsNone(val)


    @patch('builtins.input', lambda *args: 'Yes')
    def test_ask(self):
        msf_cfg.theme = msf_cfg.Theme()

        resp = msf_cfg.ask("Question", msf_cfg.check_bool)
        self.assertTrue(resp)


    @patch('builtins.input', lambda *args: '1, 2, 3')
    def test_ask_value_is_sequence(self):
        msf_cfg.theme = msf_cfg.Theme()

        resp = msf_cfg.ask("Question", msf_cfg.check_positive_int, sequence=True)
        self.assertEqual(resp, [1, 2, 3])

    @patch('builtins.input', lambda *args: 'foo')
    def test_ask_bad_value(self):
        msf_cfg.theme = msf_cfg.Theme()

        with self.catch_io(out=True):
            with self.assertRaises(RuntimeError) as ctx:
                msf_cfg.ask("Question", msf_cfg.check_bool)
        self.assertEqual(str(ctx.exception),
                         f'{msf_cfg.theme.ERROR}Too many error. Exiting{msf_cfg.theme.RESET}')

    @patch('builtins.input', lambda *args: '')
    def test_ask_use_default(self):
        msf_cfg.theme = msf_cfg.Theme()

        val = msf_cfg.ask("Question", msf_cfg.check_bool, "Yes", expected=["Yes", "No"], explanation="bla bla")
        self.assertTrue(val)


    @patch('builtins.input', lambda *args: '')
    def test_ask_default_is_sequence(self):
        msf_cfg.theme = msf_cfg.Theme()

        resp = msf_cfg.ask("Question", msf_cfg.check_positive_int, [1, 2, 3], sequence=True)
        self.assertEqual(resp, [1, 2, 3])


    @patch('builtins.input', lambda *args: 'Yes')
    def test_set_section(self):
        msf_cfg.theme = msf_cfg.Theme()
        cp = msf_cfg.ConfigParserWithComments()
        sec_name = "new_section"
        options = {'hmmer': {"question": "that is the question",
                             "validator": msf_cfg.check_choice,
                             "default": "Yes",
                             "expected": ["Yes", "no"],
                             "explanation": "an explanation"}
                   }
        defaults = MacsyDefaults()

        with self.catch_io(out=True) as out:
            msf_cfg.set_section(sec_name, options, cp, defaults, use_defaults=False)
            stdout = sys.stdout.getvalue().strip()
        self.assertTrue(cp.has_section(sec_name))
        self.assertTrue(cp.has_option(sec_name, 'hmmer'))

        self.assertEqual(stdout,
                         f"{msf_cfg.theme.SECTION}Configuring new_section options:{msf_cfg.theme.RESET}")


    @patch('builtins.input', lambda *args: 'Yes')
    def test_set_section_use_defaults(self):
        msf_cfg.theme = msf_cfg.Theme()
        cp = msf_cfg.ConfigParserWithComments()
        sec_name = "new_section"
        options = {'hmmer': {"question": "that is the question",
                             "validator": msf_cfg.check_str,
                             "default": "Yes",
                             "explanation": ""}
                   }
        defaults = MacsyDefaults(hmmer="Yes")

        with self.catch_io(out=True):
            msf_cfg.set_section(sec_name, options, cp, defaults, use_defaults=True)
            stdout = sys.stdout.getvalue().strip()
        self.assertTrue(cp.has_section(sec_name))
        self.assertFalse(cp.has_option(sec_name, 'hmmer'))

        self.assertEqual(stdout,
                         f"{msf_cfg.theme.SECTION}Configuring new_section options:{msf_cfg.theme.RESET}")

    def test_set_path_options(self):
        msf_cfg.theme = msf_cfg.Theme()
        cp = msf_cfg.ConfigParserWithComments()
        ask_ori = msf_cfg.ask
        defaults = MacsyDefaults()
        resp = ["Yes",                 # enter section ?
                None,                  # system_models_dir
                os.getcwd(),           # res_search_dir
                "res_search_suffix",   # res_search_suffix
                "res_extract_suffix",  # res_extract_suffix
                "profile_suffix"]      # profile_suffix
        g = (r for r in resp)

        def fake_ask(*args, **kwargs):
            return next(g)

        try:
            msf_cfg.ask = fake_ask
            with self.catch_io(out=True):
                msf_cfg.set_path_options(cp, defaults)
                stdout = sys.stdout.getvalue().strip()
            self.assertEqual(stdout,
                             f"{msf_cfg.theme.SECTION}Configuring directories options:{msf_cfg.theme.RESET}"
                            )
            self.assertTrue(cp.has_section("directories"))
            self.assertEqual(cp.get("directories", "res_search_dir"), os.getcwd())
            self.assertEqual(cp.get("directories", "res_search_suffix"), "res_search_suffix")
            self.assertEqual(cp.get("directories", "res_extract_suffix"), "res_extract_suffix")
            self.assertEqual(cp.get("directories", "profile_suffix"), "profile_suffix")
        finally:
            msf_cfg.ask = ask_ori


    def test_set_path_options_with_system_models_dir(self):
        msf_cfg.theme = msf_cfg.Theme()
        cp = msf_cfg.ConfigParserWithComments()
        ask_ori = msf_cfg.ask
        defaults = MacsyDefaults()

        test_dir = os.path.dirname(__file__)
        parent = os.path.join(os.path.normpath(os.path.join(test_dir, '..')))

        resp = ["Yes",                 # enter section ?
                [parent, test_dir],    # system_models_dir
                os.getcwd(),           # res_search_dir
                "res_search_suffix",   # res_search_suffix
                "res_extract_suffix",  # res_extract_suffix
                "profile_suffix"]      # profile_suffix
        g = (r for r in resp)

        def fake_ask(*args, **kwargs):
            return next(g)

        try:
            msf_cfg.ask = fake_ask
            with self.catch_io(out=True):
                msf_cfg.set_path_options(cp, defaults)
                stdout = sys.stdout.getvalue().strip()
            self.assertEqual(stdout,
                             f"{msf_cfg.theme.SECTION}Configuring directories options:{msf_cfg.theme.RESET}"
                            )
            self.assertTrue(cp.has_section("directories"))
            self.maxDiff = None
            self.assertEqual(cp.get("directories", "system_models_dir"), f"{parent}, {test_dir}")
            self.assertEqual(cp.get("directories", "res_search_dir"), os.getcwd())
            self.assertEqual(cp.get("directories", "res_search_suffix"), "res_search_suffix")
            self.assertEqual(cp.get("directories", "res_extract_suffix"), "res_extract_suffix")
            self.assertEqual(cp.get("directories", "profile_suffix"), "profile_suffix")
        finally:
            msf_cfg.ask = ask_ori


    def test_set_hmmer_options(self):
        msf_cfg.theme = msf_cfg.Theme()
        cp = msf_cfg.ConfigParserWithComments()
        ask_ori = msf_cfg.ask
        defaults = MacsyDefaults()
        resp = ["Yes",           # enter section ?
                defaults.hmmer,  # hmmer exe
                False,           # cut_ga the fake_ask do not perform casting
                0.002,           # e_value_search,
                0.003,           # i_evalue_sel,
                0.004            # coverage_profile
         ]
        g = (r for r in resp)

        def fake_ask(*args, **kwargs):
            return next(g)

        try:
            msf_cfg.ask = fake_ask
            with self.catch_io(out=True):
                msf_cfg.set_hmmer_options(cp, defaults)
                stdout = sys.stdout.getvalue().strip()
            self.assertEqual(stdout,
                             f"{msf_cfg.theme.SECTION}Configuring hmmer options:{msf_cfg.theme.RESET}")
            self.assertTrue(cp.has_section("hmmer"))
            self.assertFalse(cp.has_option("hmmer", "hmmer"))
            # all values are casted in str before inserting in ConfigParser
            self.assertEqual(cp.get("hmmer", "cut_ga"), 'False')
            self.assertEqual(cp.get("hmmer", "e_value_search"), '0.002')
            self.assertEqual(cp.get("hmmer", "i_evalue_sel"), '0.003')
            self.assertEqual(cp.get("hmmer", "coverage_profile"), '0.004')
        finally:
            msf_cfg.ask = ask_ori


    def test_set_general_options(self):
        msf_cfg.theme = msf_cfg.Theme()
        cp = msf_cfg.ConfigParserWithComments()
        ask_ori = msf_cfg.ask
        defaults = MacsyDefaults()
        resp = ["Yes",      # enter section ?
                'warning',  # log_level
                0,          # worker
                True        # mute
                ]
        g = (r for r in resp)

        def fake_ask(*args, **kwargs):
            return next(g)

        try:
            msf_cfg.ask = fake_ask
            with self.catch_io(out=True):
                msf_cfg.set_general_options(cp, defaults)
                stdout = sys.stdout.getvalue().strip()
            self.assertEqual(stdout,
                             f"{msf_cfg.theme.SECTION}Configuring general options:{msf_cfg.theme.RESET}")
            self.assertTrue(cp.has_section("general"))
            self.assertEqual(cp.get("general", "log_level"), "warning")
            self.assertEqual(cp.get("general", "worker"), "0")
            self.assertEqual(cp.get("general", "mute"), "True")
        finally:
            msf_cfg.ask = ask_ori


    def test_set_score_options(self):
        msf_cfg.theme = msf_cfg.Theme()
        cp = msf_cfg.ConfigParserWithComments()
        ask_ori = msf_cfg.ask
        defaults = MacsyDefaults()
        resp = ["Yes",  # enter section ?
                0.1,  # mandatory_weight
                0.2,  # accessory_weight
                0.3,  # exchangeable_weight
                0.4,  # redundancy_penalty
                0.5   # out_of_cluster_weight
                ]
        g = (r for r in resp)

        def fake_ask(*args, **kwargs):
            return next(g)

        try:
            msf_cfg.ask = fake_ask
            with self.catch_io(out=True):
                msf_cfg.set_score_options(cp, defaults)
                stdout = sys.stdout.getvalue().strip()
            self.assertEqual(stdout,
                             f"{msf_cfg.theme.SECTION}Configuring score_opt options:{msf_cfg.theme.RESET}")
            self.assertTrue(cp.has_section("score_opt"))
            self.assertEqual(cp.get("score_opt", "mandatory_weight"), "0.1")
            self.assertEqual(cp.get("score_opt", "accessory_weight"), "0.2")
            self.assertEqual(cp.get("score_opt", "exchangeable_weight"), "0.3")
            self.assertEqual(cp.get("score_opt", "redundancy_penalty"), "0.4")
            self.assertEqual(cp.get("score_opt", "out_of_cluster_weight"), "0.5")
        finally:
            msf_cfg.ask = ask_ori


    def test_set_base_options(self):
        msf_cfg.theme = msf_cfg.Theme()
        cp = msf_cfg.ConfigParserWithComments()
        ask_ori = msf_cfg.ask
        defaults = MacsyDefaults()
        resp = ["Yes",  # enter section ?
                "ordered",  # db_type
                "linear",  # replicon_topology
                "None",  # sequence_db
                ]
        g = (r for r in resp)

        def fake_ask(*args, **kwargs):
            return next(g)

        try:
            msf_cfg.ask = fake_ask
            with self.catch_io(out=True):
                msf_cfg.set_base_options(cp, defaults)
                stdout = sys.stdout.getvalue().strip()
            self.assertEqual(stdout,
                             f"{msf_cfg.theme.SECTION}Configuring base options:{msf_cfg.theme.RESET}")
            self.assertTrue(cp.has_section("base"))
            self.assertEqual(cp.get("base", "db_type"), "ordered")
            self.assertEqual(cp.get("base", "replicon_topology"), "linear")
            self.assertEqual(cp.get("base", "sequence_db"), "None")
        finally:
            msf_cfg.ask = ask_ori


    def test_parse_args(self):
        args = "macsyconfig --defaults"
        parsed_args = msf_cfg.parse_args(args.split()[1:])
        self.assertFalse(parsed_args.no_color)
        self.assertFalse(parsed_args.white_bg)
        self.assertTrue(parsed_args.dark_bg)
        self.assertTrue(parsed_args.defaults)

        args = "macsyconfig"
        parsed_args = msf_cfg.parse_args(args.split()[1:])
        self.assertFalse(parsed_args.no_color)
        self.assertFalse(parsed_args.white_bg)
        self.assertTrue(parsed_args.dark_bg)
        self.assertFalse(parsed_args.defaults)

        args = "macsyconfig --no-color"
        parsed_args = msf_cfg.parse_args(args.split()[1:])
        self.assertTrue(parsed_args.no_color)
        self.assertFalse(parsed_args.white_bg)
        self.assertTrue(parsed_args.dark_bg)
        self.assertFalse(parsed_args.defaults)

        args = "macsyconfig --white-bg"
        parsed_args = msf_cfg.parse_args(args.split()[1:])
        self.assertFalse(parsed_args.no_color)
        self.assertTrue(parsed_args.white_bg)
        self.assertTrue(parsed_args.dark_bg)
        self.assertFalse(parsed_args.defaults)

        real_sys_exit = argparse._sys.exit
        argparse._sys.exit = self.fake_exit
        args = "macsyconfig --white-bg --dark-bg"

        try:
            with self.catch_io(err=True):
                with self.assertRaises(TypeError) as ctx:
                    msf_cfg.parse_args(args.split()[1:])
                stderr = sys.stderr.getvalue().strip()
            self.assertEqual(str(ctx.exception), "2")  # test the return code
            self.assertEqual(stderr,
                             """usage: run_tests.py [-h] [--no-color | --white-bg | --dark-bg] [--defaults]
run_tests.py: error: argument --dark-bg: not allowed with argument --white-bg""")

        finally:
            argparse._sys.exit = real_sys_exit


    def test_functional_dark_theme(self):
        cur_dir = os.getcwd()
        tmpdir = os.path.join(tempfile.gettempdir(), 'tmp-macsyconfig')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        cf_default = MacsyDefaults(system_models_dir='')
        cf_args = argparse.Namespace(cfg_file="macsyfinder.conf")
        try:
            with self.catch_io(out=True):
                macsyconfig_args = "macsyconfig --defaults"
                msf_cfg.main(macsyconfig_args.split()[1:])
                stdout = sys.stdout.getvalue()

            expected_stdout = f"""Welcome to the MacSyFinder {msf_vers} configuration utility.

Please enter values for the following settings (just press Enter to
accept a default value, if one is given in brackets).

Configuring directories options:

Configuring hmmer options:

Configuring score_opt options:

Configuring general options:

Configuring base options:

A configuration file 'macsyfinder.conf' has been generated..
Place it in canonical location
 * in /etc/macsyfinder for system wide configuration (must named macsyfinder.conf)
 * in <VIRTUALENV>/etc if you use a virtualenv (must named macsyfinder.conf)
 * in ~/.macsyfinder for user wide configuration (must named macsyfinder.conf)
 * where you run the analysis for local configuration (must named macsyfinder.conf)
 * you can also put anywhere on the filesystems and use MACSY_CONF environment variable
   to indicate where to find it or specify it on the macsyfinder command line with option --cfg-file
   can be named as you want.


"""
            # I don't why but ansi color escape sequence are removed when I catch the stdout ???

            self.maxDiff = None
            self.assertEqual(stdout,
                             expected_stdout)
            # reparse the generated config file and check the result
            cfg = Config(cf_default, cf_args)
            for opt_name in cf_default.keys():
                if opt_name in ('inter_gene_max_space',
                                'max_nb_genes',
                                'min_genes_required',
                                'min_mandatory_genes_required',
                                'multi_loci',
                                'models_dir',
                                'out_dir'):
                    # not set in macsyconfig
                    continue
                opt_value = getattr(cfg, opt_name)()
                if opt_name == 'cfg_file':
                    self.assertEqual(opt_value, "macsyfinder.conf")
                else:
                    self.assertEqual(opt_value, cf_default[opt_name], msg=f"{opt_name}: {opt_value} {cf_default[opt_name]}")
            self.assertDictEqual(cfg.hit_weights(),
                                 {'mandatory': cf_default['mandatory_weight'],
                                  'accessory': cf_default['accessory_weight'],
                                  'neutral': cf_default['neutral_weight'],
                                  'itself': cf_default['itself_weight'],
                                  'exchangeable': cf_default['exchangeable_weight'],
                                  'out_of_cluster': cf_default['out_of_cluster_weight']
                                  })
        finally:
            os.chdir(cur_dir)
            shutil.rmtree(tmpdir)


    def test_functional_file_already_exists(self):
        cur_dir = os.getcwd()
        tmpdir = os.path.join(tempfile.gettempdir(), 'tmp-macsyconfig')
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        ask_ori = msf_cfg.ask
        real_exit = sys.exit
        sys.exit = self.fake_exit
        open('macsyfinder.conf', 'w').close()
        conf_template = self.find_data("conf_files", "macsyfinder_default.conf")
        with open(conf_template) as conf_file:
            lines = conf_file.readlines()
        # we nned to hak this line
        # as the hmmsearch path may differ on each installation
        lines[39] = f"# hmmer = {shutil.which('hmmsearch')}\n"
        expected_conf = 'macsyfinder_default.conf'
        with open(expected_conf, 'w') as conf_file:
            conf_file.write(''.join(lines))

        try:
            with self.catch_io(out=True):
                # test Abort
                macsyconfig_args = "macsyconfig --defaults"
                msf_cfg.ask = lambda *args, **kwargs: "a"
                with self.assertRaises(TypeError) as ctx:
                    msf_cfg.main(macsyconfig_args.split()[1:])
                self.assertEqual(str(ctx.exception), '1')  # check the returncode
                self.assertEqual(os.path.getsize('macsyfinder.conf'), 0)

                # test Overwrite
                macsyconfig_args = "macsyconfig --defaults"
                msf_cfg.ask = lambda *args, **kwargs: "O"
                msf_cfg.main(macsyconfig_args.split()[1:])
                # check that the produced file is equal as file in data
                self.maxDiff = None
                self.assertFileEqual(expected_conf,
                                     "macsyfinder.conf",
                                     skip_line="# system_models_dir"  # default for models dir differ on each install
                                     )
        finally:
            os.chdir(cur_dir)
            shutil.rmtree(tmpdir)
            sys.exit = real_exit
            msf_cfg.ask = ask_ori
