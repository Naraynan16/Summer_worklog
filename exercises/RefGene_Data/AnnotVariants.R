
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


partial <- as.data.frame(merge(x = partialLex,y = partialFex))

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

partialGR <- GRanges(seqnames = partial$chrom,
                     ranges = IRanges(start = partial$exStarts,
                                      end = partial$exEnds),
                     geneName = partial$geneName,
                     exonNo = partial$exonNo)
partialGR$Index <- as.character(seq_len(length(partialGR)))


intragenicGR <- GRanges(seqnames = intragenic$chrom,
                        ranges = IRanges(start = intragenic$exStarts,
                                         end = intragenic$exEnds),
                        geneName = intragenic$geneName,
                        exonNo = intragenic$exonNo)
intragenicGR$Index <- as.character(seq_len(length(intragenicGR)))

allRegions <- GRangesList(firstExonGR,
                          lastExonGR,
                          partialGR,
                          intragenicGR,
                          fiveutrGR,
                          threeutrGR)

pos <- as.vector(c("firstExon",
                   "lastExon",
                   "partial",
                   "intragenic",
                   "fiveutr",
                   "threeutr"))
#function to findOverlaps and annotate

AnnotVariants <- function(queryGR,RefGR){
  
  #check whether the queryGR has index and index it if not indexed.
  if (!"Index" %in% colnames(mcols(queryGR))) {
    queryGR$Index <- as.character(seq_len(length(queryGR)))
  }
  for (i in 1:6) {
    
    #findOverlap b/w the query and ref
    qrOverlaps <- findOverlaps(queryGR,RefGR[[i]],ignore.strand = TRUE)
    
    #store ref indexes
    ref_hits_index <- DataFrame(InRegion = RefGR[[i]][subjectHits(qrOverlaps)]$Index)
    
    #object to write the overlap regions in the query
    InRegion <- rep(NA, length(queryGR))
    
    if (nrow(ref_hits_index) > 0) {
      
      #rewrite the InRegion which contain ref indexes to geneName for annotation
      ref_hits_index$InRegion <- RefGR[[i]]$geneName[as.integer(ref_hits_index$InRegion)]
      
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
    colnames(mcols(queryGR)) <- c(oldCols,pos[i])
  }
  
  #annotate only the type
  regions <- as.data.frame(mcols(queryGR))
  type <- as.data.frame(matrix(nrow = length(rownames(regions))))
  for (i in 1:length(queryGR))
  { 
    a <- which(!is.na(regions[i,]))
    alen <- length(a)
    if(alen > 0 ){
      
      if(alen == 7){
        type[i,] <- rbind("wholeGene")
      }
      else if((!is.na(regions[i,2])) && (is.na(regions[i,c(3:7)]))){
        type[i,] <- rbind("firstExon")
      }
      else if((!is.na(regions[i,3])) && (is.na(regions[i,c(2,4:7)]))){
        type[i,] <- rbind("lastExon")
      }
      else if((!is.na(regions[i,6])) && (is.na(regions[i,c(2:5,7)]))){
        type[i,] <- rbind("fiveUtr")
      }
      else if((!is.na(regions[i,7])) && (is.na(regions[i,c(2:6)]))){
        type[i,] <- rbind("threeUtr")
      }
      else if ((!is.na(regions[i,c(2,4:5)])) && (is.na(regions[i,c(3,6:7)]))){
        type[i,]  <- rbind("partialFirstExon")
      }
      else if ((!is.na(regions[i,c(3,4:5)]) && is.na(regions[i,c(4,6:7)]))){
        type[i,]  <- rbind("partialLastExon")
      }
      else if((!is.na(regions[i,5])) && (is.na(regions[i,c(2,3,4,6,7)]))){
        type[i,] <- rbind("intragenic")
      }
      else if((!is.na(regions[i,c(2,4,5,6)])) && (is.na(regions[i,c(3,7)]))){
        type[i,] <- rbind("partialFirstexon+fiveUTR")
      }
      else if((!is.na(regions[i,c(2,4,5,7)])) && (is.na(regions[i,c(3,6)]))){
        type[i,] <- rbind("partialFirstexon+threeUTR")
      }
      else if((!is.na(regions[i,c(3,4,5,6)])) && (is.na(regions[i,c(2,7)]))){
        type[i,] <- rbind("partialLastexon+fiveUTR")
      }
      else if((!is.na(regions[i,c(3,4,5,6)])) && (is.na(regions[i,c(2,7)]))){
        type[i,] <- rbind("partialLastexon+threeUTR")
      }
      else if((!is.na(regions[i,c(2:6)])) && (is.na(regions[i,7]))){
        type[i,] <- rbind("allExons+fiveUTR")
      }
      else if((!is.na(regions[i,c(2:5,7)])) && (is.na(regions[i,6]))){
        type[i,] <- rbind("allExons+threeUTR")
      }
      else{
        type[i,] <- rbind(cbind("NoOverlaps"))
      }
    }}
  
  queryGR$type <- as.vector(unlist(type))
  return(queryGR)
}






