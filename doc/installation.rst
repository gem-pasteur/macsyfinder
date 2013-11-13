.. _installation:


************
Installation
************


TXSScan dependencies
====================
TXSScan has two dependencies, as it requires running the formatdb or the makeblastdb tools {NCBI,  #8230; Camacho, 2009 #8231} and the program Hmmer {Eddy, 2011 #6335; Eddy, 1998 #694}. Python version 2.7 is required to run TXSScan. 


Installation procedure
======================
First desarchive the source codes package, and enter the extracted folder::

  tar -xzvf txsscan-20131111.tar.gz
  cd txsscan-20131111
  
TXSScan installation follows classical "pythonic" installation procedure (see http://docs.python.org/2/install/)::

  python setup.py build
  (sudo) python setup.py install 

It is **highly recommanded** to run tests before performing the full installation::

  python setup.py build
  python setup.py -vv test 
  (sudo) python setup.py install 
  
.. note::
  ``sudo`` (*i.e.*, super-user privileges) is necesserary if you want to install it in the general file architecture.

"Classical" useful parameters for installation are supported, like for instance the --suffix option that allows to settle a different installation path (useful if you do not have super-user privileges)::

  python setup.py build --prefix /usr/local/home/bob/my_programs

will install TXSScan and required data (profiles folder and systems definition folders) in the Home directory of "bob", in the "my_progams" folder. 
  

