#!/bin/sh
./gatk-4.1.8.1/gatk AnnotateIntervals -L chr20XY.interval_list -R Homo_sapiens_assembly38.fasta -imr OVERLAPPING_ONLY -O chr20XY.annotated.tsv