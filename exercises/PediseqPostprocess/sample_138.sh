#!/bin/sh 
./gatk --java-options "-Xmx32G" PostprocessGermlineCNVCalls  --model-shard-path PediseqCNVcaller/Pediseq_1_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_2_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_3_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_4_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_5_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_6_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_7_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_8_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_9_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_10_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_11_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_12_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_13_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_14_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_15_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_16_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_17_of_18-model  --model-shard-path PediseqCNVcaller/Pediseq_18_of_18-model  --calls-shard-path PediseqCNVcaller/Pediseq_1_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_2_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_3_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_4_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_5_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_6_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_7_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_8_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_9_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_10_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_11_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_12_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_13_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_14_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_15_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_16_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_17_of_18-calls  --calls-shard-path PediseqCNVcaller/Pediseq_18_of_18-calls  --allosomal-contig X  --allosomal-contig Y  --contig-ploidy-calls ploidy-calls  --sample-index 138 --output-genotyped-intervals finalCalls/Pseq_batch16_P-Pseq-0056-P-A_1-genotyped-intervals.vcf.gz --output-genotyped-segments finalCalls/Pseq_batch16_P-Pseq-0056-P-A_1-genotyped-segments.vcf.gz  --output-denoised-copy-ratios finalCalls/Pseq_batch16_P-Pseq-0056-P-A_1-denoised-copy-ratios.tsv  --sequence-dictionary Homo_sapiens_assembly19.dict