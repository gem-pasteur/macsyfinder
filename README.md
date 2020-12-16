# MacSyFinder

[![Build Status](https://travis-ci.org/gem-pasteur/macsyfinder.svg?branch=master)](https://travis-ci.org/gem-pasteur/macsyfinder)
[![Coverage Status](https://coveralls.io/repos/github/gem-pasteur/macsyfinder/badge.svg?branch=master)](https://coveralls.io/github/gem-pasteur/macsyfinder?branch=master)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/macsyfinder)](https://pypi.org/project/macsyfinder/)
[![Open Source License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)
[![Doc](https://readthedocs.org/projects/macsyfinder/badge/?version=latest)](http://macsyfinder.readthedocs.org/en/latest/#)
[![PyPI](https://img.shields.io/pypi/v/macsyfinder)](https://pypi.org/project/macsyfinder/)
[![Docker Image Version (tag latest semver)](https://img.shields.io/docker/v/gempasteur/macsyfinder/latest)](https://hub.docker.com/repository/docker/gempasteur/macsyfinder)
![Conda](https://img.shields.io/conda/pn/bioconda/macsyfinder)

MacSyFinder - Detection of macromolecular systems in protein datasets using systems modelling and similarity search.



## Citation
Abby SS, Néron B, Ménager H, Touchon M, Rocha EPC (2014). MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems. PLoS ONE 9(10): e110726. doi:10.1371/journal.pone.0110726
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0110726


## MacSyFinder is available on pypi
 
[pypi](https://pypi.org/project/macsyfinder/)


## Installation from distribution

We encourage to install macsyfinder in a [virtualenv](https://virtualenv.pypa.io/en/latest/)

After creating a virtualenv dedicated to macsyfinder and activating it

    python3 -m venv my_project
    cd my_project
    source bin/activate

you can install macsyfinder as describe below.
    
### from pypi

    pip3 install macsyfinder==x.x

where `x.x` is the version number

### from git repository

    git clone https://github.com/gem-pasteur/macsyfinder.git
    cd macsyfinder
    pip3 install .
    
### for developers

    git clone https://github.com/gem-pasteur/macsyfinder.git
    cd macsyfinder
    pip3 install .[dev]
 
## Unit tests 

    python3 setup.py test
    
or 
    
    python3 tests/run_tests.py -vv
    
or to run a specific test

    python3 tests/run_tests.py -vv tests/test_xxx.py
        
     
### with travis-ci

[![Build Status](https://travis-ci.org/gem-pasteur/macsyfinder.svg?branch=master)](https://travis-ci.org/gem-pasteur/macsyfinder)
[![Coverage Status](https://coveralls.io/repos/github/gem-pasteur/macsyfinder/badge.svg?branch=master)](https://coveralls.io/github/gem-pasteur/macsyfinder?branch=master)

## Documentation

You will find complete documentation for setting up your project on readthedocs

[![Doc](https://readthedocs.org/projects/macsyfinder/badge/?version=latest)](http://macsyfinder.readthedocs.org/en/latest/#)

## Licence:

MacSyFinder is developed and released under [![Open Source License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Contributing 

We encourage contributions, bug report, enhancement ... 

But before to do that, we encourage to read [the contributing guide](CONTRIBUTING.md).

## Contributors

[List of all people who participated in the macsyfinder project](CONTRIBUTORS.md).

## Note

The `setsid` binary in *utils* directory is used only for functional tests on macosx. 
The binary has been build using the [setsid-macosx](https://github.com/tzvetkoff/setsid-macosx) project.