#!/bin/tcsh

set fasta = datatest/prru_psae.001.c01.fasta
python main_detect_SS.py $fasta T2SS T4P Tad >& datatest/log_prru_psae
