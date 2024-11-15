![MacSyFinder banner](./.github/logo_macsyfinder.png "MacSyFinder")

# MacSyFinder

[![Build Status](https://github.com/gem-pasteur/macsyfinder/actions/workflows/testing.yml/badge.svg?branch=master)](https://github.com/gem-pasteur/macsyfinder/actions/workflows/testing.yml)
[![codecov](https://codecov.io/gh/gem-pasteur/macsyfinder/branch/master/graph/badge.svg?token=q31HWcV3SM)](https://codecov.io/gh/gem-pasteur/macsyfinder)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/macsyfinder)](https://pypi.org/project/macsyfinder/)
[![Open Source License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)
[![Doc](https://readthedocs.org/projects/macsyfinder/badge/?version=latest)](http://macsyfinder.readthedocs.org/en/latest/#)
[![PyPI](https://img.shields.io/pypi/v/macsyfinder)](https://pypi.org/project/macsyfinder/)
[![Docker Image Version (latest semver)](https://img.shields.io/docker/v/gempasteur/macsyfinder?label=docker&sort=semver)](https://hub.docker.com/r/gempasteur/macsyfinder)
[![Conda](https://img.shields.io/conda/vn/bioconda/macsyfinder?style=plastic)](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/macsyfinder)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/gem-pasteur/macsyfinder/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/gem-pasteur/macsyfinder)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:561bfe6698ca9e58b552b4eb4e56132cac41c6f9/)](https://archive.softwareheritage.org/swh:1:dir:561bfe6698ca9e58b552b4eb4e56132cac41c6f9;origin=https://github.com/gem-pasteur/macsyfinder;visit=swh:1:snp:1bde3cb370766b10132c4e004c7cb377979928d1;anchor=swh:1:rev:868637fce184865d8e0436338af66a2648e8f6e1)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/6010/badge)](https://bestpractices.coreinfrastructure.org/projects/6010)
[![FAIR checklist badge](https://fairsoftwarechecklist.net/badge.svg)](https://fairsoftwarechecklist.net/v0.2?f=31&a=32113&i=32321&r=133)

MacSyFinder - Detection of macromolecular systems in protein datasets using systems modelling and similarity search.



## Citations

MacSyFinder v2:
Néron, Bertrand; Denise, Rémi; Coluzzi, Charles; Touchon, Marie; Rocha, Eduardo P.C.; Abby, Sophie S.
MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes.
Peer Community Journal, Volume 3 (2023), article no. e28. doi : 10.24072/pcjournal.250.
https://peercommunityjournal.org/articles/10.24072/pcjournal.250/

MacSyFinder v1:
Abby SS, Néron B, Ménager H, Touchon M, Rocha EPC (2014).
MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems.
PLoS ONE 9(10): e110726. doi:10.1371/journal.pone.0110726
http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0110726

## What new in MacSyFinder V2.x

https://macsyfinder.readthedocs.io/en/latest/user_guide/new_v2.html

## Installation

> [!IMPORTANT]
> MacSyFinder requires hmmer >= 3.1 (http://hmmer.org/).
> You need to install hmmer by yourself (except if you install macsyfinder via *conda/mamba*).
> If you are a modeler, you will need also `git`
> The other dependencies are managed by the python package manager *pip*.

### MacSyFinder is available on pypi

[![PyPI](https://img.shields.io/pypi/v/macsyfinder)](https://pypi.org/project/macsyfinder/)

### Installation from distribution

We encourage to install macsyfinder in a [virtualenv](https://virtualenv.pypa.io/en/latest/)

After creating a virtualenv dedicated to macsyfinder and activating it

    python3 -m venv my_project
    cd my_project
    source bin/activate

you can install macsyfinder as described below:

#### from pypi

    python3 -m pip install macsyfinder==x.x

where `x.x` is the version number

#### from conda/mamba

    mamba install -c bioconda macsyfinder=x.x

where `x.x` is the version number

#### from git repository

    git clone https://github.com/gem-pasteur/macsyfinder.git
    cd macsyfinder
    python3 -m pip install .

#### for modelers

https://macsyfinder.readthedocs.io/en/latest/modeler_guide/installation.html

#### for developers

https://macsyfinder.readthedocs.io/en/latest/developer_guide/installation.html

## Unit tests

    python3 setup.py test

or

    python3 tests/run_tests.py -vv

or to run a specific test

    python3 tests/run_tests.py -vv tests/test_xxx.py


### with github actions / coverage / codecov

[![Build Status](https://github.com/gem-pasteur/macsyfinder/actions/workflows/testing.yml/badge.svg?branch=master)](https://github.com/gem-pasteur/macsyfinder/actions/workflows/testing.yml)
[![codecov](https://codecov.io/gh/gem-pasteur/macsyfinder/branch/master/graph/badge.svg?token=q31HWcV3SM)](https://codecov.io/gh/gem-pasteur/macsyfinder)

## Models installation

Models are no longer shipped along macsyfinder package. To install Models you can use `macsydata`.
*macsydata* allow to manage models stored in [macsy-models](https://github.com/macsy-models).
Below some most useful commands.

  * available: List Models available on macsy-models.
  * search: Discover new packages.
  * install: Install or upgarde packages.
  * uninstall: Uninstall packages.
  * cite: How to cite a package.
  * ...

For complete documentation see
[macsydata section on readthedoc](https://macsyfinder.readthedocs.io/en/latest/user_guide/installation.html#models-installation-with-macsydata)

For models not stored in macsy-models the commands *available*, *search*, *installation from remote* or *upgrade from remote*
are **NOT** available.

For models **Not** stored in *macsy-models*, you have to manage them semi-manually.
Download the archive (do not unarchive it), then use *macsydata* for the installation.

## Documentation

You will find complete documentation for setting up your project on readthedocs

[![Doc](https://readthedocs.org/projects/macsyfinder/badge/?version=latest)](http://macsyfinder.readthedocs.org/en/latest/#)

## Example data sets

Two example datasets with command lines and expected output files are available [here](https://doi.org/10.6084/m9.figshare.21581280.v1)
and [here (for a more thorough one)](https://doi.org/10.6084/m9.figshare.21716426.v1). The 1st dataset is also
described [in the Documentation](https://macsyfinder.readthedocs.io/en/latest/user_guide/quickstart.html#an-example-data-set).

## Docker

MacSyFinder is also available as [Docker container](https://hub.docker.com/r/gempasteur/macsyfinder)

### How to use macsyfinder container with docker

The computations are performed under `msf` user in `/home/msf` inside the container.
So You have to mount a directory from the host in the container to exchange data (inputs data, and results)
from the host and the container.
The shared directory must be writable by the `msf` user or overwrite the user in the container by your id (see example below)

Furthermore the models are no longer packaged along macsyfinder. So you have to install them by yourself.
For that we provide a command line tool *macsydata* which is inspired by *pip*

    macsydata search PACKNAME
    macsydata install PACKNAME== or >=, or ... VERSION

To work with Docker you have to install models in a directory which will be mounted in the image at run time

    mkdir shared_dir
    cd shared_dir
    # install desired models in my_models
    docker run -v ${PWD}/:/home/msf -u $(id -u ${USER}):$(id -g ${USER})  gempasteur/macsyfinder:<tag> macsydata install --target /home/msf/my_models MODELS
    # run msf with these models
    docker run -v ${PWD}/:/home/msf -u $(id -u ${USER}):$(id -g ${USER})  gempasteur/macsyfinder:<tag> --db-type gembase --models-dir=/home/msf/my_models/ --models  TFF-SF Archaeal-T4P ComM MSH T2SS T4bP T4P Tad --sequence-db my_genome.fasta -w 12


### How to use with apptainer (formely Singularity)

As the *docker* image is registered in docker hub you can also use it directly with [apptainer](https://apptainer.org/docs/user/main/).
Unlike *docker* you have not to worry about shared directory, your `home` and `/tmp` are automatically shared.

    apptainer run -H ${HOME} docker://gempasteur/macsyfinder:<tag> macsydata install --target my_models MODELS
    apptainer run -H ${HOME} docker://gempasteur/macsyfinder:<tag> macsyfinder --db-type gembase --models-dir=my_models --models TFF-SF Archaeal-T4P ComM MSH T2SS T4bP T4P Tad --sequence-db my_genome.fasta -w 12

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
