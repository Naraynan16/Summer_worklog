#!/bin/sh
./gatk-4.1.8.1/gatk DetermineGermlineContigPloidy --model cohort-23wgs-20190213-contig-ploidy-model -I cvg/NA19017.tsv -O . --output-prefix ploidy-case --verbosity DEBUG