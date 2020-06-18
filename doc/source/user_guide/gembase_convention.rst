.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _gembase_convention:

**************
Gembase format
**************


In order to allow the users to run MacSyFinder on **several genomes at once**,
we propose to adopt the following convention to fulfill the requirements for the "gembase db_type".

It consists in providing for each protein, both the replicon name and a protein identifier separated by
a "_" in the first field of fasta headers. "_" are accepted in the replicon name, but not in the protein identifier.
Hence, the last "_" is the separator between the replicon name and the protein identifier.
As such, MacSyFinder will be able to treat each replicon separately to assess macromolecular systems' presence. 

For instance::

  >PlasmidA_0001 YP_003225072.1 | putative stcE protein 
  MKLKYLSCMILASLAMGAFAATAADNNSAIYFNTTQPVNDLQGGLAAEVK
  FAQSQILSAHPKEGESQQHLTSLRKSLLLVRLVKADDKTPVQVEARDAND
  KILGTLTLSPPSSLPDTVYHLDGVPADGIDFTPQNGTKKIINTVAEVNKL
  SDASGSSIKSYLANNALVEIQTANGRWIRDMYLPQGAELEGKMVRFVSYA
  GYNSTVFYGDRKVTLSVGNTLLFKYVNGQWFRSGELENNRIAYAQHTWSA
  ELPAHWIVPGLNLVIKQGNLSGSLNDINVGAPGELLLHTIDIGMLTTPRG
  RFDFAKDKEAHREYFQTIPVSRMIVNNYAPLHLKEVMLPTGTLLTDADPG
  >PlasmidA_0002 YP_003225073.1 | type II secretion protein EtpC
  MLFFLSSRRDRNLFIKDIALKMLTPNWVLCVILLIAGYQLVSVIRHFWLT
  PATSASDLSHVSVSETAVTDEHTEENFVFTLFGTASPPLSEGKVQKTTSS
  LSDDLLSGGDLDVRGILYSSVTEHSVAIFAHNNRQFSLGIGEKVPGYDAT
  ISAIFSDHIVINYQGKNASLPLRYDNPAKRNAQDDNNLIVGPVTTQANFR
  VKNIFDIMSLSPVTVNNTLSGYRLSPGKASSLFYNAGLHDNDLAVLLNGS
  ELRDTRQAKQIMKQLTELKEIKITVERDGQLYDAFIAVGEN
  ....
  >ChromosomeA_0001 YP_003573410.1 | adhesin-like protein
  MKKLFLFAALLMTGFAFYSCEDVVDNPAQDPAQSWNYSVSVKFADFDFNG
  AVDENSVPYTYKAPTTLYVLNEENTLMGTITTDAAPAIGDYGTYAGTLTG
  SIGNNLIITTKIGNDLTKQDGTLKSAIENGIVQTAEVPIKIYNANSGTLT
  TASAKMDNTAAIAYTSLGYIKGGDKILFVEGNQTFEWTVNEEFDPYTSTD
  LYIALPMNTDPETEYTISSDSKDGYTRGGTFKLADYPTLAAGKVSNYIGG
  IPFIQTGVDLTKWDAYMRTDPNNTWYMNNINNGWPATFSQEVEDGKSFIV
  TQSGPTLDSLNVVVGGVTGKEVNVTLNNIRLGKDRSINIGDKHGWVEYDG
  THDIYGWGAKANVTLIGENECETLYIQCPATKKGEGTLNYKNLSIDSYGS
  >ChromosomeA_0020 YP_003573411.1 | hypothetical protein
  MKRIVLITLVSILTTFQAIAQVANGFYRVQNNASSRYITLRDNAVGTVDY
  SSTNVDLSNIVTWSGFDKVKSNPASIIYVEQHDSKYDLKVQGTGIYAITG
  GRTYLELRPKDSGYILAVTYNGMEGRLYDSEEDVDGEGYVKRSGNSAYQY
  WSFIPVDTENNYIGLQPTVQVGDNYYGTLYASYPFKAASSGIKFYYVDAI
  ....
  >NC_001548_0015 YP_003225080.1 | type II secretion protein EtpJ  (translation)
  MSQQRVKGFTLLEMLLALAVFAALSISAFQVLQSGIRAHELSQDKVRRLA
  ELQRGGSQIERDLMQMIPRHSRGSEGLLLAAPHLLKSDDWGISFTRNSWL
  NPAGMLPRPELQWVGYRLRQQKLERLSYFYVDHPSGIAPDVRVVLEGVHA
  FRLRFFVNGTWQARWDSTSILPQAVEVTLVMDDFAELTRLFLVSKETAE

This input file contains 3 replicons: PlasmidA (which 2 first protein identifiers are 0001 and 0002),
ChromosomeA (which 2 first protein identifiers are 0001 and 0020) and NC_001548 (which first protein identifier is 0015).
MacSyFinder search results will thus be reported for each of these three replicons. 

.. _topology-files:

**************
Topology files
**************

To be able to attribute a topology per replicon/genome when using the Gembase format,
we propose the user to build a "topology file" in the form of a tabular file with two columns separated by a ":".
The 1st column is the replicon name, and the 2nd the corresponding topology. Comments can be written after a "#".

For example::

  # comment line
  PlasmidA : circular
  ChromosomeA : linear
  ChromosomeB : circular
  
.. note::
    A topology file can be specified on the command-line with the ``--topology-file`` parameter.
    
