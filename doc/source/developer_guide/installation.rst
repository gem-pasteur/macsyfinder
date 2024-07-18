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

==================================
MacSyFinder Installation procedure
==================================

Installation steps:
=======================

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
