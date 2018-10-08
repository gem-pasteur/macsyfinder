```
data_set_3
    base            => Contains Genome and Index.
    models          => Contains Models definition and Profiles.
    results         => Contains output produced by MacSyFinder when processing this dataset.
                       Command line details:
                       ./bin/macsyfinder --out-dir=tests/data/data_set_3/results \
                                         --relative-path \
                                         --sequence-db=tests/data/data_set_3/base/VICH001.B.00001.C001.prt \
                                         --db-type=gembase \
                                         --models-dir=tests/data/data_set_3/models \
                                         --models set_1 T2SS T4P
```
