# Secretion Systems (TXSScan)

This folder contains the definitions for 15 types of secretion systems or bacterial appendages (T1SS, T2SS, T3SS, T4P, pT4SSt, pT4SSi, T5aSS, T5bSS, T5bSS, T6SSi, T6SSii, T6SSiii, Flagellum, Tad, T9SS).

The basic command is:

    macsyfinder secretion_system \
                --db-type ordered_replicon \
                -d Macsyfinder_models/Data/TXSS/DEF \
                -p Macsyfinder_models/Data/TXSS/HMM \
                --profile-suffix .hmm \
                --sequence-db my_protein_file

where `secretion_system` is one of the list above, or `all` (see [MacSyFinder's documentation](http://macsyfinder.readthedocs.io/en/latest/))

# Download

To download only this module: [click here](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/gem-pasteur/Macsyfinder_models/tree/master/Data/TXSS)

# References

- Abby Sophie S., Néron Bertrand, Ménager Hervé, Touchon Marie, Rocha Eduardo P. C. (2014). MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems. In PLoS ONE, 9 (10), pp. e110726. [doi:10.1371/journal.pone.0110726](http://dx.doi.org/10.1371/journal.pone.0110726)

- Abby Sophie S., Cury Jean, Guglielmini Julien, Néron Bertrand, Touchon Marie, Rocha Eduardo P. C. (2016). Identification of protein secretion systems in bacterial genomes. In Scientific Reports, 6, pp. 23080. [doi:10.1038/srep23080](http://dx.doi.org/10.1038/srep23080)
