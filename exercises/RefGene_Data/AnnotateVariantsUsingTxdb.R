dat <-as.data.frame(read.csv(file = "test.exomeDepth.rawCalls.txt",header = T,sep ="\t"))
dat$chromosome <- paste0("chr",dat$chromosome)

library(GenomicRanges)
datGr <- GRanges(seqnames = dat$chromosome,
                 ranges = IRanges(start = dat$start,end = dat$end))
library(VariantAnnotation)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
allVaraints_txdb <- locateVariants(query = datGr,subject = txdb_hg19,AllVariants())
                
