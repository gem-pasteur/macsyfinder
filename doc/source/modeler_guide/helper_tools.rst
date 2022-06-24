.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2022  Institut Pasteur (Paris),and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _helper_tool:

***********
Helper Tool
***********

.. _macsyprofile:

macsyprofile
============

To help develop new models we provide the tool `macsyprofile` which is to be used as post treatement.

It is ran over a previous macsyfinder analysis:

   * it extracts from raw HMMER output files the hits and computes the profile coverage for each of them.
   * it enables to filter the hits in a user-defined manner, to test other values of filtering parameters than those used with the MacSyFinder run.
   * it writes down the results in a file in `tsv` format `hmm_coverage.tsv`.

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
      --index-dir INDEX_DIR
                            Specifies the path to a directory to store/read the sequence index
                            when the sequence-db dir is not writable.
      -f, --force           force to write output even the file already exists
                            (overwrite it).
      -V, --version         show program's version number and exit
      -v, --verbosity       Increases the verbosity level. There are 4 levels:
                            Error messages (default), Warning (-v), Info (-vv) and
                            Debug.(-vvv)
      --mute                Mute the log on stdout. (continue to log on
                            macsyfinder.log) (default: False)

    For more details, visit the MacSyFinder website and see the MacSyFinder documentation.

For instance:

.. code-block:: shell

    >macsyprofile  macsyfinder-2021XXXX_XX-XX-XX

will analyse the HMMER raw outputs stored in `macsyfinder-2021XXXX_XX-XX-XX/hmmer_results` directory
and the results will be stored in `macsyfinder-2021XXXX_XX-XX-XX/hmm_coverage.tsv` file


Setting filtering parameters
----------------------------

This helper tool is designed to help the user test the relevance of the HMM profiles used, what filtering parameters for HMMER to be used, and understand why some components might be unexpectedly missing from the MacSyFinder results.
This can thus help to improve the models - for instance for the genomic location parameters (is a component not found cause it should be listed as a `loner`?).

Therefore by default, the filtering parameters are very loose so that most hits found with HMMER will be reported, even the weakest ones.

However, it is possible to filter hits to be extracted based on the profile coverage with `--coverage-profile` or the i-evalue (`--i-evalue-sel`) to be a bit more stringent.

Also, it is possible to use the `--best-hits` in order to report only the best hit for a given protein sequence when several profiles were matching hit.


Using patterns with "--pattern"
-------------------------------

If in `previous_run/hmmer_results` you have the following files:

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


But you are interested only in ComM family genes, you can specify the option ``--pattern 'ComM*'``
For instance:

.. code-block:: text

    >macsyprofile --pattern 'ComM*'  macsyfinder-2021XXXX_XX-XX-XX
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comB.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comC.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comEA.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comEB.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comEC.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comGA.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comGB.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comGC.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comGD.search_hmm.out
    parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comGE.search_hmm.out
    found 79 hits
    result is in 'macsyfinder-2021XXXX_XX-XX-XX/hmm_coverage.tsv'

.. note::

    The patterns available are the `glob` patterns (the jokers usable with unix `ls` command )

    .. code-block:: text

        >macsyprofile --pattern 'ComM_com?C' -f macsyfinder-2021XXXX_XX-XX-XX
        parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comEC.search_hmm.out
        parsing macsyfinder-2021XXXX_XX-XX-XX/hmmer_results/ComM_comGC.search_hmm.out
        found 16 hits
        result is in 'macsyfinder-2021XXXX_XX-XX-XX/hmm_coverage.tsv'


A useful example for modellers?
-------------------------------

.. code-block:: text

    >macsyprofile --best-hits i_eval --i-evalue-sel 0.001 --coverage-profile 0.5 -o msf_GCF_003149495.1_ASM314949v1_tff-sf/hmm_coverage_best-hits_ieval_default_filter_MSF.tsv msf_GCF_003149495.1_ASM314949v1_tff-sf
    found 221 hits
    result is in 'msf_GCF_003149495.1_ASM314949v1_tff-sf/hmm_coverage_best-hits_ieval_default_filter_MSF.tsv'

This command line might be useful to macsy-models modellers, as it consists in extracting all relevant hits that are used by
the MacSyFinder engine to search systems, when using the default parameters:

- the proteins are assigned with their best hits (i-evalue based) when they match several profiles (`--best-hits i_eval` option)
- the default filtering parameters (i-evalue and profile coverage) are used (`--i-evalue-sel` and `--coverage-profile` options)

By using this command line that lists all hits available for MacSyFinder to search for systems, one could be interested in
comparing this list to the list of hits that end in being assigned to systems (listed e.g. in best_solution.tsv).
This can help to determine why a component is missing from a system: is it because there are no good hits for it, or
is it because it does not comply to the co-localization rules defined in the systems' model?




Parsing macsyprofile outputs
----------------------------

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

        systems = pd.read_csv("path/to/hmm_coverage.tsv", sep='\t', comment='#')




.. warning::

    The `macsyprofile` tool is not compliant with results produced with `macsyfinder v1`.
    If you get ``Cannot find models in conf file XXX. May be these results have been generated with an old version of macsyfinder.``
    Check the configuration file, if `[models]` section contains ``models_1 = XXX YYY`` remove the `_1` from models
    ``models = XXX YYY``
