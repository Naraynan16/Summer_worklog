---
title: summer_worklog
author: sathya
date: 08/08/2020
output: html_document
---
started with learning the shing package in R and creating apps

date :26/06/2020 -
 Todat I learnt what is git ? and basic git operations
 git init / branch / add /commit / merge / remote add / push
 created a TCGA_expression repository and  submitted the files into a repository using git in command-line

 date :30/06/2020 -
 learnt how to insert a image in rmd file.
 also solved the yaml header problem in the previous files.

 date: 01/07/2020 
 learnt the refgene data and viewed them in the ucsc genome browser
 
 date:05/07/2020
 extracted the coding transcripts from the refgene
 uploaded it into the git. 

date:10/07/2020
wrote a code to extract the utr from the coding transcripts 
and resolved the git issue. reuploaded the summer worklog repo

date:11/07/2020
wrote code to extract the introns fromm the refgene

date:12/07/2020
wrote a shiny app to view all the intervals of particular gene

date:13/07/2020
did some correction in the refgene_transcript_extract.R.
i.e replace  na.omit() with column specific omition
read about the data.table package and used fread to read a file

date:14/07/2020
revised the genomic ranges 
and started doing exercise 1- annotation of exomdepth data with refgene.

date:15/07/2020
Started with find the overlap regions of the exomeDepth data
using GenomicRanges against the genomic produced earlier using refgene data
Faced issue after annotation.

date:16/07/2020
used a reference code to understand the annotation of the Exomedepth data
https://github.com/drramki-chop/edm/blob/master/R/variant.annotations.R
wrote a code to annotate data using Txdb.

date:17-21/07/2020
learnt ways to annotate the data and wrote a code to annotate variants
in the exome data.
Modified the code to get precise info about the variants.

date:01-08/08/2020
started learning about the gcnv pipeline by GATK
https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants
used tutorial data to learn and produce data.
Formed a frame-work on the thesis project

