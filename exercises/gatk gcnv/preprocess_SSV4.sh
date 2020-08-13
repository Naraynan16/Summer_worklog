#!/bin/sh
./gatk-4.1.8.1/gatk PreprocessIntervals -R Homo_sapiens_assembly38.fasta --padding 0 -L SSV4.bed -imr OVERLAPPING_ONLY -O preProcessed_SSV4.interval_list