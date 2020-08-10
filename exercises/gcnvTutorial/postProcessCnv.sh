#!/bin/sh
./gatk-4.1.8.1/gatk PostprocessGermlineCNVCalls --model-shard-path cohort24-twelve/cohort24-twelve_1of2-model --model-shard-path cohort24-twelve/cohort24-twelve_2of2-model --calls-shard-path cohort24-twelve/cohort24-twelve_1of2-calls --calls-shard-path cohort24-twelve/cohort24-twelve_2of2-calls --allosomal-contig chrX --allosomal-contig chrY --contig-ploidy-calls ploidy-calls --sample-index 19 --output-genotyped-intervals genotyped-intervals-cohort24-twelve-NA19017.vcf.gz --output-genotyped-segments genotyped-segments-cohort24-twelve-NA19017.vcf.gz --output-denoised-copy-ratios denoised_copy_ratios-cohort24-twelve-NA19017.tsv --sequence-dictionary Homo_sapiens_assembly38.dict