#!/bin/tcsh

set bin=$TXSSCAN_HOME/bin/txsscan
$bin ../test/datatest/prru_psae.001.c01.fasta T2SS -d ../data/DEF -p ../data/profiles -r ./datatest/res_search -v
