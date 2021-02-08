.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020  Institut Pasteur (Paris),and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _helper_tool:

***********
Helper Tool
***********

macsyprofile
============

To help to develop new model we provide a tool `macsyprofile` which is used as post treatement.
It is run on previous macsyfinder analysis, it extract from raw hmmer files the hits and compute the profile coverage
for each of them.
Then it write down the results in a file in `tsv` format.

.. code-block:: text

    usage: macsyprofile [-h] [--coverage-profile COVERAGE_PROFILE]
                        [--i-evalue-sel I_EVALUE_SEL]
                        [--best-hits {score,i_eval,profile_coverage}] [-p PATTERN]
                        [-o OUT] [-f] [-V] [-v] [--mute]
                        previous_run

         *            *               *                   * *       *
    *           *               *   *   *  *    **                *
      **     *    *   *  *     *                    *               *
         __  __    *       ____ *      ____        **   __ _ _  *
        |  \/  | __ _  ___/ ___| _   _|  _ \ _ __ ___  / _(_) | ___
        | |\/| |/ _` |/ __\___ \| | | | |_) | '__/ _ \| |_| | |/ _ \
        | |  | | (_| | (__ ___) | |_| |  __/| | | (_) |  _| | |  __/
        |_|  |_|\__,_|\___|____/ \__, |_|   |_|  \___/|_| |_|_|\___|
                *                |___/    *                   *
     *      *   * *     *   **         *   *  *           *
      *      *         *        *    *              *
                 *                           *  *           *     *

    MacSyProfile - MacSyFinder profile helper tool

    positional arguments:
      previous_run          The path to a macsyfinder results directory.

    optional arguments:
      -h, --help            show this help message and exit
      --coverage-profile COVERAGE_PROFILE
                            Minimal profile coverage required for the hit
                            alignment with the profile to allow the hit selection
                            for systems detection. (default no threshold)
      --i-evalue-sel I_EVALUE_SEL
                            Maximal independent e-value for Hmmer hits to be
                            selected for systems detection. (default: no selection
                            based on i-evalue)
      --best-hits {score,i_eval,profile_coverage}
                            If several hits match the same replicon, same gene.
                            Select only the best one (based on best 'score' or
                            'i_evalue' or 'profile_coverage')
      -p PATTERN, --pattern PATTERN
                            pattern to filter the hmm files to analyse.
      -o OUT, --out OUT     the path to a file to write results.
      -f, --force           force to write output even the file already exists
                            (overwrite it).
      -V, --version         show program's version number and exit
      -v, --verbosity       Increases the verbosity level. There are 4 levels:
                            Error messages (default), Warning (-v), Info (-vv) and
                            Debug.(-vvv)
      --mute                Mute the log on stdout. (continue to log on
                            macsyfinder.log) (default: False)

    For more details, visit the MacSyFinder website and see the MacSyFinder documentation.



--pattern example
-----------------

If in previous_run hmmer_results you have the following files

.. code-block:: text

    previous_run/hmmer_results/Archaeal-T4P_arCOG11238.search_hmm.out
    previous_run/hmmer_results/Archaeal-T4P_arCOG11520.search_hmm.out
    previous_run/hmmer_results/Archaeal-T4P_arCOG11777.search_hmm.out
    previous_run/hmmer_results/Archaeal-T4P_arCOG11778.search_hmm.out
    previous_run/hmmer_results/Archaeal-T4P_arCOG11936.search_hmm.out
    previous_run/hmmer_results/Archaeal-T4P_arCOG14515.search_hmm.out
    previous_run/hmmer_results/ComM_comC.search_hmm.out
    previous_run/hmmer_results/ComM_comEB.search_hmm.out
    previous_run/hmmer_results/ComM_comEC.search_hmm.out
    previous_run/hmmer_results/ComM_comGA.search_hmm.out
    previous_run/hmmer_results/ComM_comGB.search_hmm.out
    previous_run/hmmer_results/ComM_comGC.search_hmm.out
    previous_run/hmmer_results/ComM_comGD.search_hmm.out
    previous_run/hmmer_results/ComM_comGE.search_hmm.out
    previous_run/hmmer_results/MSH_mshA.search_hmm.out
    previous_run/hmmer_results/MSH_mshB.search_hmm.out
    previous_run/hmmer_results/MSH_mshC.search_hmm.out


But you are interested only in ComM family genes, you can specify the option --pattern 'ComM*'.

.. code-block:: text

    previous_run/hmmer_results/ComM_comC.search_hmm.out
    previous_run/hmmer_results/ComM_comEB.search_hmm.out
    previous_run/hmmer_results/ComM_comEC.search_hmm.out
    previous_run/hmmer_results/ComM_comGA.search_hmm.out
    previous_run/hmmer_results/ComM_comGB.search_hmm.out
    previous_run/hmmer_results/ComM_comGC.search_hmm.out
    previous_run/hmmer_results/ComM_comGD.search_hmm.out
    previous_run/hmmer_results/ComM_comGE.search_hmm.out

The patterns availables are the `glob` patterns (the jokers usable with unix `ls` command )

The `macsyprofile` output is a tabulated separated values (`.tsv`) files
The first lines which are comments (starting with '#') display the tool version
and the complete command line used. Then follow the results.
The first line of results is a header line.

.. code-block:: text

    # macsyprofile 2.0rc1
    # macsyprofile --pattern ComM* --coverage-profile 0.5 macsyfinder-20201202_15-17-46/
    hit_id  replicon_name   position_hit    hit_sequence_length     gene_name       i_eval  score   profile_coverage        sequence_coverage       begin   end
    GCF_000006745_021980    GCF_000006745   2198    291     ComM_comC       2.500e-40       136.400 0.942   0.708   62      267
    GCF_000006745_007650    GCF_000006745   765     253     ComM_comC       9.600e-31       105.100 0.937   0.798   43      244
    ...

.. note::
    This file can be easily parsed using the Python `pandas <https://pandas.pydata.org/>`_ library. ::

        import pandas as pd

        systems = pd.read_cvs("path/to/hmm_coverage.tsv", sep='\t', comment='#')

