.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2024 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _dev_installation:


************
Installation
************

MacSyFinder works with models for macromolecular systems that are not shipped with it,
you have to install them separately. See the :ref:`macsydata section <modeler_macsydata>` below.
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

Then run `macsyfinder` or `macsydata`.


.. note::

    from 2.1.4 version, *MacSyFinder* has adopted `ruff <https://docs.astral.sh/ruff/>`_ as linter
    and *pre-commit* to ensure the coding style.
    please read `CONTRIBUTING.md <https://github.com/gem-pasteur/macsyfinder/blob/master/CONTRIBUTING.md>`_ guide lines.


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
    test_average_wholeness (tests.test_solution.SolutionTest.test_average_wholeness) ... ok
    test_gt (tests.test_solution.SolutionTest.test_gt) ... ok
    test_hits_number (tests.test_solution.SolutionTest.test_hits_number) ... ok
    test_hits_positions (tests.test_solution.SolutionTest.test_hits_positions) ... ok
    test_iteration (tests.test_solution.SolutionTest.test_iteration) ... ok
    test_lt (tests.test_solution.SolutionTest.test_lt) ... ok
    test_score (tests.test_solution.SolutionTest.test_score) ... ok
    test_sorting (tests.test_solution.SolutionTest.test_sorting) ... ok
    test_systems (tests.test_solution.SolutionTest.test_systems) ... ok
    test_get_def_to_detect (tests.test_utils.TestUtils.test_get_def_to_detect) ... ok
    test_get_replicon_names_bad_type (tests.test_utils.TestUtils.test_get_replicon_names_bad_type) ... ok
    test_get_replicon_names_gembase (tests.test_utils.TestUtils.test_get_replicon_names_gembase) ... ok
    test_get_replicon_names_ordered (tests.test_utils.TestUtils.test_get_replicon_names_ordered) ... ok
    test_get_replicon_names_unordered (tests.test_utils.TestUtils.test_get_replicon_names_unordered) ... ok
    test_parse_time (tests.test_utils.TestUtils.test_parse_time) ... ok
    test_threads_available (tests.test_utils.TestUtils.test_threads_available) ... ok

    ----------------------------------------------------------------------
    Ran 548 tests in 34.265s

    OK

The tests must be in python file (`.py`) starting with with `test\_` \
It's possible to specify one or several test files, one module, or one class in a module or a method in a Test class.

Test the `test_package` module

.. code-block:: shell

    python -m unittest -vv tests.test_package

.. code-block:: text

    test_init (tests.test_package.TestLocalModelIndex.test_init) ... ok
    test_repos_url (tests.test_package.TestLocalModelIndex.test_repos_url) ... ok
    test_check (tests.test_package.TestPackage.test_check) ... ok
    test_check_bad_metadata (tests.test_package.TestPackage.test_check_bad_metadata) ... ok
    test_check_dir_in_profile (tests.test_package.TestPackage.test_check_dir_in_profile) ... ok

    ...

    test_list_package_vers (tests.test_package.TestRemoteModelIndex.test_list_package_vers) ... ok
    test_list_packages (tests.test_package.TestRemoteModelIndex.test_list_packages) ... ok
    test_remote_exists (tests.test_package.TestRemoteModelIndex.test_remote_exists) ... ok
    test_repos_url (tests.test_package.TestRemoteModelIndex.test_repos_url) ... ok
    test_unarchive (tests.test_package.TestRemoteModelIndex.test_unarchive) ... ok
    test_url_json (tests.test_package.TestRemoteModelIndex.test_url_json) ... ok
    test_url_json_reach_limit (tests.test_package.TestRemoteModelIndex.test_url_json_reach_limit) ... ok

    ----------------------------------------------------------------------
    Ran 56 tests in 0.242s

    OK

Test only the class `TestPackage` (this module contains 3 classes)

.. code-block:: shell

    python -m unittest -vv tests.test_package.TestPackage

.. code-block:: text

    test_check (tests.test_package.TestPackage.test_check) ... ok
    test_check_bad_metadata (tests.test_package.TestPackage.test_check_bad_metadata) ... ok
    test_check_dir_in_profile (tests.test_package.TestPackage.test_check_dir_in_profile) ... ok
    test_check_empty_profile (tests.test_package.TestPackage.test_check_empty_profile) ... ok
    test_check_metadata (tests.test_package.TestPackage.test_check_metadata) ... ok
    test_check_metadata_no_cite (tests.test_package.TestPackage.test_check_metadata_no_cite) ... ok

    ...

    test_metadata (tests.test_package.TestPackage.test_metadata) ... ok
    test_profile_with_bad_ext (tests.test_package.TestPackage.test_profile_with_bad_ext) ... ok

    ----------------------------------------------------------------------
    Ran 42 tests in 0.151s

    OK

Test only the method `test_metadata` from the test Class `TestPackage` in module `test_package`

.. code-block:: shell

    python -m unittest -vv tests.test_package.TestPackage.test_metadata

.. code-block:: text

    test_metadata (tests.test_package.TestPackage.test_metadata) ... ok

    ----------------------------------------------------------------------
    Ran 1 test in 0.005s

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

    test_lt (tests.test_solution.SolutionTest.test_lt) ... ok
    test_score (tests.test_solution.SolutionTest.test_score) ... ok
    test_sorting (tests.test_solution.SolutionTest.test_sorting) ... ok
    test_systems (tests.test_solution.SolutionTest.test_systems) ... ok
    test_get_def_to_detect (tests.test_utils.TestUtils.test_get_def_to_detect) ... ok
    test_get_replicon_names_bad_type (tests.test_utils.TestUtils.test_get_replicon_names_bad_type) ... ok
    test_get_replicon_names_gembase (tests.test_utils.TestUtils.test_get_replicon_names_gembase) ... ok
    test_get_replicon_names_ordered (tests.test_utils.TestUtils.test_get_replicon_names_ordered) ... ok
    test_get_replicon_names_unordered (tests.test_utils.TestUtils.test_get_replicon_names_unordered) ... ok
    test_parse_time (tests.test_utils.TestUtils.test_parse_time) ... ok
    test_threads_available (tests.test_utils.TestUtils.test_threads_available) ... ok

    ----------------------------------------------------------------------
    Ran 548 tests in 34.485s

    OK

Then display a report

.. code-block:: shell

    coverage report


.. code-block:: text

    Name                                     Stmts   Miss Branch BrPart  Cover
    --------------------------------------------------------------------------
    macsypy/__init__.py                         57      4     12      1    93%
    macsypy/cluster.py                         218      1    122      2    99%
    macsypy/config.py                          390     10    160      9    97%
    macsypy/database.py                        204      3     74      1    99%
    macsypy/definition_parser.py               220      3     76      2    98%
    macsypy/error.py                             9      0      0      0   100%
    macsypy/gene.py                            143      4     48      3    96%
    macsypy/hit.py                             199      1     83      2    99%
    macsypy/licenses.py                         14      1      2      1    88%
    macsypy/metadata.py                        126      0     76      2    99%
    macsypy/model.py                           127      0     58      0   100%
    macsypy/model_conf_parser.py                62      0     16      0   100%
    macsypy/package.py                         326      9    144      5    96%
    macsypy/profile.py                         116      7     34      1    95%
    macsypy/registries.py                      189      5     80      6    96%
    macsypy/report.py                          122      0     44      2    99%
    macsypy/scripts/__init__.py                  0      0      0      0   100%
    macsypy/scripts/macsy_gembase_split.py      94      2     30      3    96%
    macsypy/scripts/macsy_merge_results.py     201     15    102      8    91%
    macsypy/scripts/macsyconfig.py             242      3    108      4    98%
    macsypy/scripts/macsydata.py               690     62    209     16    90%
    macsypy/scripts/macsyfinder.py             551     23    243     10    96%
    macsypy/scripts/macsyprofile.py            247      5     86      6    97%
    macsypy/search_genes.py                     80      7     24      3    90%
    macsypy/serialization.py                   133      4     72      3    97%
    macsypy/solution.py                         97      0     74      0   100%
    macsypy/system.py                          397      3    204      0    99%
    macsypy/utils.py                            69      0     33      1    99%
    --------------------------------------------------------------------------
    TOTAL                                     5323    172   2214     91    96%

or generate a html report

.. code-block:: shell

    coverage html

.. code-block:: text

    Wrote HTML report to htmlcov/index.html

The results are in the `htmlcov` directory. With you favourite web browser, open the `index.html` file.
for more options please refer to the `coverage documentation <https://coverage.readthedocs.io/en/latest/>`_ .
