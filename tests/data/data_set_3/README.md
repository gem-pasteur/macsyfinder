```
data_set_3
    base            => Contains Genome and Index.
    models          => Contains Models definition and Profiles.
    results         => Contains output produced by MacSyFinder when processing this dataset.
                       Command line details:
                       ./bin/macsyfinder --out-dir=tests/data/data_set_2/results \
                                         --relative-path \
                                         --sequence-db=tests/data/data_set_2/base/test.fa \
                                         --db-type=gembase \
                                         --models-dir=tests/data/data_set_2/models \
                                         --models set_1 T9SS T3SS T4SS_typeI
```
