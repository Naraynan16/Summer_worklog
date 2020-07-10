refgene <- read.table("C:/Users/Narayanan/Downloads/refGene.txt")
header <- c("id", "transcript_id", "chrom", "txstrand", "txStart", "txEnd", "cdsStart",
            "cdsEnd", "exonCount", "exStarts", "exEnds", "uid", "geneName", "cdsStartStat","cdsEndStat","exonFrames")
colnames(refgene) <- header

#extract coding transcripts

g1 <- grep("NM_", refgene$transcript_id)
codingTx <- refgene[g1, ]
write.csv(file = "coding_transcripts.csv", x = codingTx,row.names = F)
allExons<- data.frame()

library(dplyr)
library(tidyr)
allExons <- codingTx %>% separate_rows(exStarts,exEnds)

write.csv(x = allExons,file = "individual_transcripts.csv",row.names = F)
