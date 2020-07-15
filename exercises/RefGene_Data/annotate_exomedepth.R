dat <-as.data.frame(read.csv(file = "test.exomeDepth.rawCalls.txt",header = T,sep ="\t"))
dat$chromosome <- paste0("chr",dat$chromosome)
introns <- read.csv(file = "refgene_introns.csv",header = T)
exons <- read.csv(file = "individual_transcripts.csv",header = T)
utr1 <- as.data.frame(read.csv(file = "utr1.csv",header = T))
utr1 <- utr1[!utr1$start== "-",]
utr2 <- as.data.frame(read.csv(file = "utr2.csv",header = T))
utr2 <- utr2[!utr2$start== "-",]

dat$intervalType <- NA
for (i in 1:dim(dat)[1]){
  chr <- dat[i,7]
  exlist <- exons %>% filter(exons$chrom == chr)
  exOverlaps <- exlist[between(x = dat[i,5],lower = exlist$exStarts,upper = exlist$exEnds,incbounds = F),c(1,6)]
  if(dim(exOverlaps)[1]>0){
    dat$intervalType <- "exonic"
  }
 inlist <- introns %>% filter(introns$chrom == chr)
 inOverlaps<-inlist[between(x = dat[i,5],lower = inlist$start,upper = inlist$end,incbounds = F),c(1,7)]
 if(dim(inOverlaps)[1]>0){
   dat$intervalType <- "intronic"
 }                             
  utlist1 <- utr1 %>% filter(utr1$chrom == chr)
  utrOverlaps1 <- utlist1[between(x = dat[i,5],lower = utlist1$start,upper = utlist1$start,incbounds = F),c(4,7,1)]
  if(dim(utrOverlaps1)[1]>0){
    dat$intervalType <- utrOverlaps1[,3]
  }                            
  utlist2 <- utr2 %>% filter(utr2$chrom == chr)
   utrOverlaps2 <- utlist2[between(x = dat[i,5],lower = utlist2$start,upper = utlist2$start,incbounds = F),c(4,7,1)]
  if(dim(utrOverlaps2)[1]>0){
    dat$intervalType <- utrOverlaps2[,4]
  }                           
}


