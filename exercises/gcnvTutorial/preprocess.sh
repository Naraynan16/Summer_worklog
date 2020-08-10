#!/bin/sh
./gatk-4.1.8.1/gatk PreprocessIntervals -R Homo_sapiens_assembly38.fasta --padding 0 -L gcnv-chr20XY-contig.list -imr OVERLAPPING_ONLY -O chr20XY.interval_list
