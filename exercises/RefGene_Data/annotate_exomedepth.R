dat <-as.data.frame(read.csv(file = "test.exomeDepth.rawCalls.txt",header = T,sep ="\t"))
dat$chromosome <- paste0("chr",dat$chromosome)
introns <- read.csv(file = "refgene_introns.csv",header = T)
exons <- read.csv(file = "individual_transcripts.csv",header = T)
utr1 <- as.data.frame(read.csv(file = "utr1.csv",header = T))
utr1 <- utr1[!utr1$start== "-",]
utr2 <- as.data.frame(read.csv(file = "utr2.csv",header = T))
utr2 <- utr2[!utr2$start== "-",]

datGr <- GRanges(seqnames = dat$chromosome,
                 ranges = IRanges(start = dat$start,end = dat$end,type= dat$type))

exonsGr <- GRanges(seqnames = exons$chrom,
                   ranges = IRanges(start = exons$exStarts,end = exons$exEnds,geneName= exons$geneName)) 

intronGr<- GRanges(seqnames = introns$chrom,
                   ranges = IRanges(start = introns$start,end = introns$end,geneName=introns$geneName))

utrGr1<- GRanges(seqnames = utr1$chrom,
                 ranges = IRanges(start = type.convert(utr1$start),
                                  end = type.convert(utr1$end),
                                  type = utr1$type))
utrGr2<- GRanges(seqnames = utr2$chrom,
                 ranges = IRanges(start = as.numeric(utr2$start),end = as.numeric(utr2$end),
                                  type = utr2$type))


exonicOverlaps <- findOverlaps(datGr,exonsGr)
exSubset <- subsetByOverlaps(exonsGr,datGr)
intronicOverlaps <- findOverlaps(datGr,intronGr)
utr1Overlaps <- findOverlaps(datGr,utrGr1)
utr2Overlaps <- findOverlaps(datGr,utrGr2)

subsetByOverlaps(exonsGr,datGr)


