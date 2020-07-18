
#input
dat <-as.data.frame(read.csv(file ="test.exomeDepth.rawCalls.txt",header = T,sep ="\t"))
dat$chromosome <- paste0("chr",dat$chromosome)

#GRanges for dat
library(GenomicRanges)

datGr <- GRanges(seqnames = dat$chromosome,
                 ranges = IRanges(start = dat$start,end = dat$end))

#index for easy element access and internal matching between GRanges
datGr$Index <- as.character(seq_len(length(datGr)))

#datas for annotation

introns <- read.csv(file = "refgene_introns.csv",header = T)
exons <- read.csv(file = "individual_transcripts.csv",header = T)
fiveutr <- as.data.frame(read.csv(file = "fiveUTR.csv",header = T))
threeutr <- as.data.frame(read.csv(file = "threeUTR.csv",header = T))

#firstExon
library(dplyr)
library(tidyr)

firstExon <- as.data.frame(exons %>%
                             group_by(transcript_id) %>%
                             filter(row_number() == 1 ))

#lastExon and  omitting genes with one exon only since they are in the firstExon object.

lastExon <-as.data.frame(exons[exons$exonCount != 1, ] %>%
                           group_by(transcript_id) %>%
                           filter(row_number()== n()))

#intragenic is variant not in the first and last exon
tmp <- as.data.frame(exons %>%
                       group_by(transcript_id) %>%
                       filter(row_number() != 1)) 

intragenic <- as.data.frame(tmp %>%
                              group_by(transcript_id) %>%
                              filter(row_number() != n())) 


#partial gene with first exon overlap 

partialFex <- as.data.frame(exons %>%
                              group_by(transcript_id) %>%
                              filter(row_number() != n()))

#partial gene with last exon overlap

partialLex <- as.data.frame(exons %>%
                              group_by(transcript_id) %>%
                              filter(row_number() != 1))

#making GRanges for the genomicRegions

firstExonGR <- GRanges(seqnames = firstExon$chrom,
                       ranges = IRanges(start = firstExon$exStarts,
                                        end = firstExon$exEnds),
                       geneName = firstExon$geneName,
                       exonNo = firstExon$exonNo)
firstExonGR$Index <- as.character(seq_len(length(firstExonGR)))

lastExonGR <- GRanges(seqnames = lastExon$chrom,
                      ranges = IRanges(start = lastExon$exStarts,
                                       end = lastExon$exEnds),
                      geneName = lastExon$geneName,
                      exonNo = lastExon$exonNo)
lastExonGR$Index <- as.character(seq_len(length(lastExonGR)))

threeutrGR <- GRanges(seqnames = threeutr$chrom,
                      ranges = IRanges(start = threeutr$start,
                                       end = threeutr$end),
                      geneName = threeutr$geneName,
                      type = threeutr$type)
threeutrGR$Index <- as.character(seq_len(length(threeutrGR)))

fiveutrGR <- GRanges(seqnames = fiveutr$chrom,
                     ranges = IRanges(start = fiveutr$start,
                                      end = fiveutr$end),
                     geneName = fiveutr$geneName,
                     type = fiveutr$type)
fiveutrGR$Index <- as.character(seq_len(length(fiveutrGR)))

partialFexGR <- GRanges(seqnames = partialFex$chrom,
                        ranges = IRanges(start = partialFex$exStarts,
                                         end = partialFex$exEnds),
                        geneName = partialFex$geneName,
                        exonNo = partialFex$exonNo)
partialFexGR$Index <- as.character(seq_len(length(partialFexGR)))

partialLexGR <- GRanges(seqnames = partialLex$chrom,
                        ranges = IRanges(start = partialLex$exStarts,
                                         end = partialLex$exEnds),
                        geneName = partialLex$geneName,
                        exonNo = partialLex$exonNo)
partialLexGR$Index <- as.character(seq_len(length(partialLexGR)))

intragenicGR <- GRanges(seqnames = intragenic$chrom,
                        ranges = IRanges(start = intragenic$exStarts,
                                         end = intragenic$exEnds),
                        geneName = intragenic$geneName,
                        exonNo = intragenic$exonNo)
intragenicGR$Index <- as.character(seq_len(length(intragenicGR)))

intronsGR <- GRanges(seqnames = introns$chrom,
                     ranges = IRanges(start = introns$start,
                                      end = introns$end),
                     geneName = introns$geneName,
                     exonNo = introns$intron_no)
intronGr$Index <- as.character(seq_len(length(intronGr)))

#function to findOverlaps and annotate

AnnotVariants <- function(queryGR,RefGR,pos){
  
  #check whether the queryGR has index and index it if not indexed.
  if (!"Index" %in% colnames(mcols(queryGR))) {
    queryGR$Index <- as.character(seq_len(length(queryGR)))
  }
      #findOverlap b/w the query and ref
  qrOverlaps <- findOverlaps(queryGR,RefGR,ignore.strand = TRUE)
  
  #store ref indexes
  ref_hits_index <- DataFrame(InRegion = RefGR[subjectHits(qrOverlaps)]$Index)
  
  #object to write the overlap regions in the query
  InRegion <- rep(NA, length(queryGR))
  
  if (nrow(ref_hits_index) > 0) {
    
    #rewrite the InRegion which contain ref indexes to geneName for annotation
    ref_hits_index$InRegion <- RefGR$geneName[as.integer(ref_hits_index$InRegion)]
    
    #func to collapse genes for a particular hit  
    clps <- function(x) paste0(x, collapse = "|")
    
    #combine hits obtained for each query
    
    ref_hits_index_combined <- aggregate(
      ref_hits_index,
      list(Index = as.character(queryGR$Index[queryHits(qrOverlaps)])),
      FUN =clps)
    
    #write the obtained annotation to an object to add it to the query
    InRegion[base::match(ref_hits_index_combined$Index, 
                         queryGR$Index)] <- ref_hits_index_combined$InRegion
  }
  
  #add annotation to query and 
  oldCols<- colnames(mcols(queryGR))
  queryGR$newcolumn <- InRegion
  colnames(mcols(queryGR)) <- c(oldCols,pos)
  return(queryGR)
}

datGr <- AnnotVariants(queryGR = datGr,RefGR = firstExonGR,pos = "FirstExon")
datGr <- AnnotVariants(queryGR = datGr,RefGR = lastExonGR,pos = "LastExon")
datGr <- AnnotVariants(queryGR = datGr,RefGR = threeutrGR,pos = "threeUTR")
datGr <- AnnotVariants(queryGR = datGr,RefGR = fiveutrGR,pos = "fiveUTR")
datGr <- AnnotVariants(queryGR = datGr,RefGR = intragenicGR,pos = "intragenic")
datGr <- AnnotVariants(queryGR = datGr,RefGR = partialFexGR,pos = "partial_with_firstExon")
datGr <- AnnotVariants(queryGR = datGr,RefGR = partialLexGR,pos = "partial_with_lastExon")
datGr <- AnnotVariants(queryGR = datGr,RefGR = intronGr,pos = "intron")



