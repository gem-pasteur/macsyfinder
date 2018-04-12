INSTALLATION INSTRUCTIONS FOR MacSyFinder and MacSyView
=======================================================


1 - Requirements :
==================

for MacSyFinder :
-----------------
- Any machine running a unix-like Operating System.
- Python, (>=2.7 & <3.0)
- makeblastdb (>= 2.2.28) or formatdb (>=2.2.26)
   + makeblastdb is a distribute with blastdb 
     (http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
   + formatdb 
- hmmsearch 3.1 (http://hmmer.janelia.org/)


for MacSyView :
---------------
- chromium or firefox web browser.

2 - Archive overview:
=====================

* bin => contain the executables (macsyfinder)
* data => the hmm profiles and the systems definitions
* doc => the documentation in html and pdf
* etc => a template of macsyfinder configuration file
* test => all needed for unit tests
* macsypy => the macsyfinder python library
* setup.py => the installation script

3 - Installation steps:
=======================

3.1 - Make sure every required dependence/software is present.
--------------------------------------------------------------

By default MacSyFinder will search if `makeblastdb` is on your PATH if not it will looking for `formatdb`.
if neither `makeblastdb` nor `formatdb` are in your path you have to set the absolute path toward one of 
these tools in the configuration file.
by default MacSyFinder will try to use `hmmsearch` in your PATH, if `hmmsearch` is not in the PATH, 
you have to set the absolute path to the `hmmsearch` in configuration.
If the tools are not in the path, some test will be skipped and a warning will be raised.


3.2 - Perform the installation.
-------------------------------

###Uncompress the archive


```bash
tar -xzf macsyfinder-x.x.tar.gz
cd macsyfinder-x.x
```

###Build MacSyFinder

macsyfinder installation follows classical "pythonic" installation procedure (see http://docs.python.org/2/install/):

```bash
python setup.py build
```

###Test MacSyFinder

To test the libraries you just build:

```bash
python setup.py test -vv
```

###Install

```bash
sudo python setup.py install
```

If you have not the privileges to perform a system wide installation,
you can use a [virtual environment](https://virtualenv.pypa.io/en/stable/).

```bash
vitualenv macsyfinder
source macsyfinder/bin/activate
mkdir macsyfinder/src
cd macsyfinder/src
```

and follow the "regular" procedure:
* uncompress the archive
* build MacSyFinder
* test MacSyFinder

for installation in virtualenv don't use `--prefix` nor `sudo`

```bash
python setup.py install
```

To exit the virtualenv just execute the `deactivate` command.
To run `macsyfinder`, you need to activate the virtualenv

```bash
source macsyfinder/bin/activate
```

Then run macsyfinder/macsyview


You can access all general cmd with  
```bash
python setup.py --help
```
or specific options with 

```bash
python setup.py cmd --help
```


4 - Configuration.
==================
macsyfinder have some options that can be set in configuration file 
 + system wide $prefix/etc/macsyfinder/macsyfinder.conf
 + user wide ~/.macsyfinder/macsyfinder.conf
 + for project ./macsyfinder.conf
or on the command line.

for all options see documentation 
 + ($prefix or /usr)/local/share/doc/macsyfinder/html/index.html
 + ($prefix or /usr)/local/share/doc/macsyfinder/pdf/Macsyfinder.pdf
 
 
5 - Uninstalation.
==================

To uninstall MacSyFinder (the last version installed), run:
(sudo) python setup.py uninstall


To get started after installation, see the Running MacSyFinder section in the
 User's Guide (Userguide.pdf)