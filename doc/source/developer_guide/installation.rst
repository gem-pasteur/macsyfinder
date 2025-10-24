.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2025 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details
    MacSyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _dev_installation:


************
Installation
************

MacSyFinder works with models for macromolecular systems that are not shipped with it,
you have to install them separately. See the :ref:`msf_data section <modeler_msf_data>` below.
We also provide container so you can use macsyfinder directly.

.. dev_dependencies:

========================
MacSyFinder dependencies
========================
**Python version >=3.10** is required to run MacSyFinder: https://docs.python.org/3.10/index.html

MacSyFinder has one program dependency:

 - the *Hmmer* program, version 3.1 or greater (http://hmmer.org/).

The *hmmsearch* program should be installed (*e.g.*, in the PATH) in order to use MacSyFinder.
Otherwise, the paths to this executable must be specified in the command-line:
see the :ref:`command-line options <hmmer-options>`.


MacSyFinder also relies on some Python library dependencies:

 - MacSyLib
 - colorlog
 - colorama
 - pyyaml
 - packaging
 - networkx
 - pandas
 - GitPython
 - sphinx
 - sphinx_rtd_theme
 - sphinx-autodoc-typehints
 - sphinxcontrib-svg2pdfconverter
 - coverage
 - build
 - ruff
 - pre-commit

These dependencies will be automatically retrieved and installed when using `pip` for installation (see below).


.. dev_install:

===============================================
MacSyFinder Installation and testing procedures
===============================================

Installation steps:
===================

Make sure every required dependency/software is present.
--------------------------------------------------------

By default MacSyFinder will try to use `hmmsearch` in your PATH. If `hmmsearch` is not in the PATH,
you have to set the absolute path to `hmmsearch` in a :ref:`configuration file <config-definition-label>`
or in the :ref:`command-line <hmmer-options>` upon execution.
If the tools are not in the path, some test will be skipped and a warning will be raised.


installation in a virtualenv
""""""""""""""""""""""""""""

.. code-block:: bash

    # create a new virtaulenv
    python3 -m venv macsyfinder
    # activate it
    cd macsyfinder
    source bin/activate
    # clone/install the project in editable mode
    git clone
    cd macsyfinder
    python3 -m pip install -e .[dev]
    # install tools to ensure coding style
    pre-commit install

To exit the virtualenv just execute the `deactivate` command.

.. code-block:: bash

    source macsyfinder/bin/activate

Then run `macsyfinder` or `msf_data`.


.. note::

    from 2.1.4 version, *MacSyFinder* has adopted `ruff <https://docs.astral.sh/ruff/>`_ as linter
    and *pre-commit* to ensure the coding style.
    please read `CONTRIBUTING.md <https://github.com/gem-pasteur/macsyfinder/blob/main/CONTRIBUTING.md>`_ guide lines.


.. dev_testing:

Testing
=======

MacSyFinder project use `unittest` framework (included in the standard library) to test the code.

All tests stuff is in `tests` directory.

* The data directory contains data needed by the tests
* in the __init__.py file a MacsyTest class is defined and should be the base of all testcase use in the project
* each test_*.py represent a file containing unit or functional tests

To run all the tests (in the virtualenv)

.. code-block:: shell

    python -m unittest discover

To increase verbosity of output

.. code-block:: shell

    python -m unittest discover -vv

.. code-block:: text

    ...
    test_check_choice (tests.test_macsyconfig.TestMacsyconfig.test_check_choice) ... ok
    test_check_dir (tests.test_macsyconfig.TestMacsyconfig.test_check_dir) ... ok
    test_check_exe (tests.test_macsyconfig.TestMacsyconfig.test_check_exe) ... ok
    test_check_file (tests.test_macsyconfig.TestMacsyconfig.test_check_file) ... ok
    test_check_float (tests.test_macsyconfig.TestMacsyconfig.test_check_float) ... ok
    test_check_positive_int (tests.test_macsyconfig.TestMacsyconfig.test_check_positive_int) ... ok
    test_check_str (tests.test_macsyconfig.TestMacsyconfig.test_check_str) ... ok
    test_functional_dark_theme (tests.test_macsyconfig.TestMacsyconfig.test_functional_dark_theme) ... ok
    test_functional_file_already_exists (tests.test_macsyconfig.TestMacsyconfig.test_functional_file_already_exists) ... ok
    test_parse_args (tests.test_macsyconfig.TestMacsyconfig.test_parse_args) ... ok
    test_set_base_options (tests.test_macsyconfig.TestMacsyconfig.test_set_base_options) ... ok
    test_set_general_options (tests.test_macsyconfig.TestMacsyconfig.test_set_general_options) ... ok
    test_set_hmmer_options (tests.test_macsyconfig.TestMacsyconfig.test_set_hmmer_options) ... ok
    test_set_path_options (tests.test_macsyconfig.TestMacsyconfig.test_set_path_options) ... ok
    test_set_path_options_with_system_models_dir (tests.test_macsyconfig.TestMacsyconfig.test_set_path_options_with_system_models_dir) ... ok
    test_set_score_options (tests.test_macsyconfig.TestMacsyconfig.test_set_score_options) ... ok
    test_set_section (tests.test_macsyconfig.TestMacsyconfig.test_set_section) ... ok
    test_set_section_use_defaults (tests.test_macsyconfig.TestMacsyconfig.test_set_section_use_defaults) ... ok
    test_db_type_set_to_gembase (tests.test_macsyfinder.TestMacsyfinder.test_db_type_set_to_gembase) ... ok
    test_get_version_message (tests.test_macsyfinder.TestMacsyfinder.test_get_version_message) ... ok
    test_list_models (tests.test_macsyfinder.TestMacsyfinder.test_list_models) ... ok
    test_list_models_no_permissions (tests.test_macsyfinder.TestMacsyfinder.test_list_models_no_permissions) ... ok
    test_parse_args (tests.test_macsyfinder.TestMacsyfinder.test_parse_args) ... ok
    test_search_systems (tests.test_macsyfinder.TestMacsyfinder.test_search_systems) ... ok
    test_search_systems_unordered (tests.test_macsyfinder.TestMacsyfinder.test_search_systems_unordered) ... ok
    test_cmde (tests.test_scripts.Test_macsyfinder.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_config.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_data.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_merge.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_profile.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_split.test_cmde) ... ok

    ----------------------------------------------------------------------
    Ran 85 tests in 8.339s

    OK

The tests must be in python file (`.py`) starting with with `test\_` \
It's possible to specify one or several test files, one module, or one class in a module or a method in a Test class.

Test the `test_macsyfinder` module

.. code-block:: shell

    python -m unittest -vv tests.test_macsyfinder

.. code-block:: text

    test_db_type_set_to_gembase (tests.test_macsyfinder.TestMacsyfinder.test_db_type_set_to_gembase) ... ok
    test_get_version_message (tests.test_macsyfinder.TestMacsyfinder.test_get_version_message) ... ok
    test_list_models (tests.test_macsyfinder.TestMacsyfinder.test_list_models) ... ok
    test_list_models_no_permissions (tests.test_macsyfinder.TestMacsyfinder.test_list_models_no_permissions) ... ok
    test_parse_args (tests.test_macsyfinder.TestMacsyfinder.test_parse_args) ... ok
    test_search_systems (tests.test_macsyfinder.TestMacsyfinder.test_search_systems) ... ok
    test_search_systems_unordered (tests.test_macsyfinder.TestMacsyfinder.test_search_systems_unordered) ... ok

    ----------------------------------------------------------------------
    Ran 7 tests in 1.773s

    OK

Test only the class `TestMacsyfinder`

.. code-block:: shell

    python -m unittest -vv tests.test_macsyfinder.TestMacsyfinder

.. code-block:: text

    test_db_type_set_to_gembase (tests.test_macsyfinder.TestMacsyfinder.test_db_type_set_to_gembase) ... ok
    test_get_version_message (tests.test_macsyfinder.TestMacsyfinder.test_get_version_message) ... ok
    test_list_models (tests.test_macsyfinder.TestMacsyfinder.test_list_models) ... ok
    test_list_models_no_permissions (tests.test_macsyfinder.TestMacsyfinder.test_list_models_no_permissions) ... ok
    test_parse_args (tests.test_macsyfinder.TestMacsyfinder.test_parse_args) ... ok
    test_search_systems (tests.test_macsyfinder.TestMacsyfinder.test_search_systems) ... ok
    test_search_systems_unordered (tests.test_macsyfinder.TestMacsyfinder.test_search_systems_unordered) ... ok

    ----------------------------------------------------------------------
    Ran 7 tests in 1.862s

    OK

Test only the method `test_list_models_no_permissions` from the test Class `TestMacsyfinder` in module `test_macsyfinder`

.. code-block:: shell

    python -m unittest -vv tests.test_macsyfinder.TestMacsyfinder.test_list_models_no_permissions

.. code-block:: text

    test_list_models_no_permissions (tests.test_macsyfinder.TestMacsyfinder.test_list_models_no_permissions) ... ok

    ----------------------------------------------------------------------
    Ran 1 test in 0.000s

    OK


Coverage
========

To compute the tests coverage, we use the `coverage <https://pypi.org/project/coverage/>`_ package.
The package is automatically installed if you have installed `macsyfinder` with the `dev` target see :ref:`installation <dev_installation>`
The coverage package is setup in the `pyproject.toml` configuration file

To compute the coverage

.. code-block:: shell

    coverage run

.. code-block:: text

    ...

    test_set_score_options (tests.test_macsyconfig.TestMacsyconfig.test_set_score_options) ... ok
    test_set_section (tests.test_macsyconfig.TestMacsyconfig.test_set_section) ... ok
    test_set_section_use_defaults (tests.test_macsyconfig.TestMacsyconfig.test_set_section_use_defaults) ... ok
    test_db_type_set_to_gembase (tests.test_macsyfinder.TestMacsyfinder.test_db_type_set_to_gembase) ... ok
    test_get_version_message (tests.test_macsyfinder.TestMacsyfinder.test_get_version_message) ... ok
    test_list_models (tests.test_macsyfinder.TestMacsyfinder.test_list_models) ... ok
    test_list_models_no_permissions (tests.test_macsyfinder.TestMacsyfinder.test_list_models_no_permissions) ... ok
    test_parse_args (tests.test_macsyfinder.TestMacsyfinder.test_parse_args) ... ok
    test_search_systems (tests.test_macsyfinder.TestMacsyfinder.test_search_systems) ... ok
    test_search_systems_unordered (tests.test_macsyfinder.TestMacsyfinder.test_search_systems_unordered) ... ok
    test_cmde (tests.test_scripts.Test_macsyfinder.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_config.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_data.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_merge.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_profile.test_cmde) ... ok
    test_cmde (tests.test_scripts.Test_msf_split.test_cmde) ... ok

    ----------------------------------------------------------------------
    Ran 85 tests in 8.339s

    OK

Then display a report

.. code-block:: shell

    coverage report


.. code-block:: text

    Name                                             Stmts   Miss Branch BrPart  Cover
    ----------------------------------------------------------------------------------
    src/macsyfinder/__init__.py                         23      2      0      0    91%
    src/macsyfinder/io.py                               11      0      4      0   100%
    src/macsyfinder/scripts/macsy_gembase_split.py      95      2     22      3    96%
    src/macsyfinder/scripts/macsy_merge_results.py     203     15     72      8    90%
    src/macsyfinder/scripts/macsyconfig.py             245      5     88      4    97%
    src/macsyfinder/scripts/msf.py                     237      8     48      4    96%
    ----------------------------------------------------------------------------------
    TOTAL

or generate a html report

.. code-block:: shell

    coverage html

.. code-block:: text

    Wrote HTML report to htmlcov/index.html

The results are in the `htmlcov` directory. With you favourite web browser, open the `index.html` file.
for more options please refer to the `coverage documentation <https://coverage.readthedocs.io/en/latest/>`_ .
