#!/bin/sh
gatk CollectReadCounts -L chr20sub.interval_list -R Homo_sapiens_assembly38.fasta -imr OVERLAPPING_ONLY -I NA19017.chr20sub.bam --format TSV -O NA19017.tsv 