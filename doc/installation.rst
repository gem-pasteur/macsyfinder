.. _installation:


************
Installation
************


TXSScan dependencies
====================
TXSScan has two dependencies, as it requires running the formatdb (>=2.2.26)or 
the makeblastdb (2.2.28) tools provided with the Blast suite of programs 
(http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 
and the program Hmmer version 3.0 (http://hmmer.janelia.org/). 
Python version 2.7 is required to run TXSScan. 
Thus formatdb or makeblastdb and hmmsearch programs must be installed (*e.g.*, in the PATH) in order to use TXSScan. 


Installation procedure
======================
First unarchive the source codes package, and enter the extracted folder::

  tar -xzvf txsscan-x.x.tar.gz
  cd txsscan-x.x
  
TXSScan installation follows classical "pythonic" installation procedure (see http://docs.python.org/2/install/)::

  python setup.py build
  (sudo) python setup.py install 

It is **highly recommanded** to run tests before performing the full installation::

  python setup.py build
  python setup.py test -vv
  (sudo) python setup.py install 
  
.. note::
  super-user privileges (*i.e.*, ``sudo``) are necesserary if you want to install it in the general file architecture.
  
  
.. note::
  If you have not the privileges, or you don't want, to install TXSScan in the python libraries of your system, 
  you can install TXSScan in a virtualenv (http://www.virtualenv.org/)

Procedures specific to TXSScan can be used instead of default. Please run the command for full options::
  

  python setup.py --help

The main ones are::
 
  python setup.py install --prefix /usr/local/home/bob/my_programs # Specifies an installation path

which will install TXSScan and required data (profiles folder and systems definition folders) in the Home directory of "bob", in the "my_progams" folder (useful if you do not have super-user privileges)

To uninstall TXSScan (the last version installed), run::

  (sudo) python setup.py uninstall 

