# ```data_set_8```

## Base

* gembase.fa
* ordered_replicon.fa
* unordered.fa
* unordered_replicon.fa

## Model

'multi-loci' flag has been set to TRUE for T9SS model 
and set to FALSE for T3SS and T4SS_typeI models.

## Directory tree description

```
data_set_2
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
                         --out-dir=tests/data/data_set_2/results \
                         --relative-path \
                         --sequence-db=tests/data/data_set_2/base/test.fa \
                         --db-type=gembase \
                         --models-dir=tests/data/data_set_2/models \
                         --models set_1 T9SS T3SS T4SS_typeI
```
