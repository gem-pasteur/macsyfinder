# MacSyFinder

[![pipeline status](https://gitlab.pasteur.fr/gem/MacSyFinder/badges/master/pipeline.svg)](https://gitlab.pasteur.fr/gem/MacSyFinder/commits/master)
[![coverage report](https://gitlab.pasteur.fr/gem/MacSyFinder/badges/master/coverage.svg)](https://gitlab.pasteur.fr/gem/MacSyFinder/commits/master)
[![Doc](https://img.shields.io/badge/docs-passed-brightgreen.svg)](http://gem.pages.pasteur.fr/MacSyFinder/)

MacSyFinder - Detection of macromolecular systems in protein datasets using systems modelling and similarity search.



## Citation
Abby SS, Néron B, Ménager H, Touchon M, Rocha EPC (2014). MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems. PLoS ONE 9(10): e110726. doi:10.1371/journal.pone.0110726
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0110726


## Download distribution
 
[ ![Download](https://api.bintray.com/packages/gem-pasteur/MacSyFinder/macsyfinder/images/download.svg) ](https://bintray.com/gem-pasteur/MacSyFinder/macsyfinder/_latestVersion)


## Installation from distribution

We encourage to install macsyfinder in a [virtualenv](https://virtualenv.pypa.io/en/latest/)

After creating a virtualenv dedicated to macsyfinder and activating it

    python -m venv my_project
    cd my_project
    source bin/activate

you can install macsyfinder as describe below.
    
### from distribution

    pip install macsyfinder-xxx.tar.gz

### from git repository

    git clone https://gitlab.pasteur.fr/gem/MacSyFinder.git
    cd macsyfinder
    pip install .
    
### for developers

    git clone https://gitlab.pasteur.fr/gem/MacSyFinder.git
    cd macsyfinder
    pip install .[dev]
 
## Unit tests 

    python setup.py test
    
or 
    
    python tests/run_tests.py -vv
    
or to run a specific test

    python tests/run_tests.py -vv tests/test_xxx.py
        
     
### with gitlab-ci

[![pipeline status](https://gitlab.pasteur.fr/gem/MacSyFinder/badges/master/pipeline.svg)](https://gitlab.pasteur.fr/gem/MacSyFinder/commits/master)  
[![coverage report](https://gitlab.pasteur.fr/gem/MacSyFinder/badges/master/coverage.svg)](http://gem.pages.pasteur.fr/MacSyFinder/coverage/index.html)

## Documentation

You will find complete documentation for setting up your project.
for the development version on gitlab pages

[![Doc](https://img.shields.io/badge/docs-passed-brightgreen.svg)](http://gem.pages.pasteur.fr/MacSyFinder/)


for public version on readthedocs

[![Doc](https://readthedocs.org/projects/macsyfinder/badge/?version=latest)](http://macsyfinder.readthedocs.org/en/latest/#)

## Licence:

MacSyFinder is developed and released under [open source licence GPLv3](https://opensource.org/licenses/GPL-3.0)

## Contributing 

We encourage contributions, bug report, enhancement ... 

But before to do that, we encourage to read [the contributing guide](CONTRIBUTING.md).

## Note

The `setsid` binary in *utils* directory is used only for functional tests on macosx. 
The binary has been build using the [setsid-macosx](https://github.com/tzvetkoff/setsid-macosx) project.