# ```data_set_9```

## Base

unordered_replicon.fa

## Model

'multi-loci' flag has been set to FALSE for T9SS, T3SS
and T4SS_typeI models.

## Directory tree description

```
data_set_9
    base            => Contains Genome and Index.
        test.fa.ori => This file was used to produce the output
                       in the 'results' folder (during the run, 
                       test.fa.ori was named test.fa).
                       This file is not used during test execution.
        test.fa     => Manually modifed version of test.fa.ori file, 
                       created after the production of 'results' folder.
                       This file is used during test execution.
    models          => Contains Models definition and Profiles.
    results         => Contains output produced by MacSyFinder
                       when processing this dataset.
                       Command line details:
                       ./bin/macsyfinder \
                         --out-dir=tests/data/data_set_9/results \
                         --relative-path \
                         --sequence-db=tests/data/data_set_9/base/unordered_replicon.fa \
                         --db-type=unordered_replicon \
                         --models-dir=tests/data/data_set_9/models \
                         --models set_1 T9SS T3SS T4SS_typeI
```
