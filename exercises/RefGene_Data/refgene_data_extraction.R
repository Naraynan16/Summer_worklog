library(data.table)
refgene <- as.data.frame(fread("C:/Users/Narayanan/Downloads/refgene.txt.gz"))
colnames(refgene) <- c("id", "transcript_id", "chrom", "txstrand", "txStart", "txEnd", "cdsStart",
                      "cdsEnd", "exonCount", "exStarts", "exEnds", "uid", "geneName", "cdsStartStat",
                      "cdsEndStat","exonFrames")

#extract coding transcripts

g1 <- grep("NM_", refgene$transcript_id)
codingTx <- refgene[g1, ]
write.csv(file = "coding_transcripts.csv", x = codingTx,row.names = F)
allExons<- data.frame()

library(dplyr)
library(tidyr)

set1  <- codingTx[,c(2:4,9,10,11,13)]
allExons <- set1 %>% separate_rows(exStarts,exEnds,convert = T)
allExons <-  allExons[complete.cases(allExons$exStarts),]
allExons <- as.data.frame(allExons %>% group_by(transcript_id) %>% mutate(exonNo = paste0("exon_",row_number())))
write.csv(x = allExons,file = "individual_transcripts.csv",row.names = F)
