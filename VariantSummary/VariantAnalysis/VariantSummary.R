library(VariantAnnotation)
library(GenomicRanges)

#targeted interval list used in sequencing
SSV4_interval <-as.data.frame(read.delim(file = "SSV4.bed",header = F,sep = "\t"))
colnames(SSV4_interval)<- c("chr","start","end")
SSV4_intervalRanges <- GRanges(seqnames = SSV4_interval$chr,
                               ranges = IRanges(start = SSV4_interval$start,
                                                end = SSV4_interval$end))
#sampleinfo
sampleName <- as.data.frame(read.csv(file = "samples.txt",header = F))
fileName <- paste0("C:/Users/Narayanan/Desktop/summer_worklog/PediseqVcfs/",sampleName$V1,"-genotyped-segments.vcf")

#per sample details
sampleDetails <- as.data.frame(matrix(nrow = dim(sampleName)[1],ncol = 16))
colnames(sampleDetails) <- c("SampleName","nCNV","meanCNVSize","minCNVsize","maxCNVsize",
                             "nDel","meanDelSize","Del.Std.deviation","minDelSize","maxDelSize",
                             "nDup","meanDupSize","Dup.Std.deviation","minDupSize","maxDupSize",
                             "MeanTargetExons")

#combine all cnvs for all samples
allCNVs <- data.frame()

for (i in 1:length(fileName)){
 #read vcf
  vcf <- readVcf(file = fileName[i],genome = "hg19")

  #extract data from vcf
  vcfRanges <- rowRanges(vcf)
  type <- unlist(alt(vcf))
  id <- names(vcfRanges)
  
  #a data frame is created to arrange CNV data
  CNV <- data.frame(matrix(unlist(strsplit(id,"_")), nrow=length(id), byrow=T))
  CNV <- CNV[,-1]
  colnames(CNV) <- c("chr","start","end")
  CNV$start <- as.numeric(CNV$start)
  CNV$end <- as.numeric(CNV$end)
  CNV$type <- type
  CNV <- CNV[!CNV$type=="", ]
  nCNVs <- as.numeric(dim(CNV)[1])
  CNV$sampleName <- rep(sampleName[i, 1],dim(CNV)[1])
  cnvRanges <- GRanges(seqnames = CNV$chr,
                       ranges = IRanges(start = CNV$start,end = CNV$end),type=CNV$type)
  
  CNV$size <- as.numeric(width(cnvRanges))
  
  #number of dels and dups
  nDels <- as.numeric(dim(CNV[CNV$type=="<DEL>", ])[1])
  nDups <-as.numeric(dim(CNV[CNV$type=="<DUP>", ])[1])
  
  #find overlaps to get number of target exons
  CNV$nTargetExons <- as.numeric(countOverlaps(query =cnvRanges ,subject = SSV4_intervalRanges))
    
  sampleDetails[i,] <-rbind(c(sampleName[i,1],
                              nCNVs,
                              mean(CNV$size),
                              min(CNV$size),
                              max(CNV$size),
                              nDels,
                              mean(CNV[CNV$type=="<DEL>",6]),
                              sd(CNV[CNV$type=="<DEL>",6]),
                              min(CNV[CNV$type=="<DEL>",6]),
                              max(CNV[CNV$type=="<DEL>",6]),
                              nDups,
                              mean(CNV[CNV$type=="<DUP>",6]),
                              sd(CNV[CNV$type=="<DUP>",6]),
                              min(CNV[CNV$type=="<DUP>",6]),
                              max(CNV[CNV$type=="<DUP>",6]),
                              mean(CNV$nTargetExons)))
 
   allCNVs <- rbind(allCNVs,CNV)
  
}

allDels <-allCNVs[allCNVs$type=="<DEL>", ]
allDups <-allCNVs[allCNVs$type=="<DUP>", ]


sampleDetails$meanCNVSize <- as.numeric(sampleDetails$meanCNVSize)
sampleDetails$minCNVsize <- as.numeric(sampleDetails$minCNVsize)
sampleDetails$maxCNVsize <- as.numeric(sampleDetails$maxCNVsize)
sampleDetails$meanDelSize <- as.numeric(sampleDetails$meanDelSize)
sampleDetails$Del.Std.deviation <- as.numeric(sampleDetails$Del.Std.deviation)
sampleDetails$minDelSize <- as.numeric(sampleDetails$minDelSize)
sampleDetails$maxDelSize <- as.numeric(sampleDetails$maxDelSize)
sampleDetails$meanDupSize <- as.numeric(sampleDetails$meanDupSize)
sampleDetails$Dup.Std.deviation <- as.numeric(sampleDetails$Dup.Std.deviation)
sampleDetails$minDupSize <- as.numeric(sampleDetails$minDupSize)
sampleDetails$maxDupSize <- as.numeric(sampleDetails$maxDupSize)
sampleDetails$MeanTargetExons <- as.numeric(sampleDetails$MeanTargetExons)



#Summary of allCNVs

totalCNV <- dim(allCNVs)[1]
mean_nCNVs <-  mean(as.numeric(sampleDetails$nCNV)) 
std_devCNV <- sd(as.numeric(sampleDetails$nCNV))

Samples_below_two_std_dev <- sampleDetails[sampleDetails$nCNV < mean_nCNVs-(2*std_devCNV) , ]
Samples_above_two_std_dev <- sampleDetails[sampleDetails$nCNV > mean_nCNVs+(2*std_devCNV) , ]

meanCNVsize <- mean(allCNVs$size)
minCNVSize <- min(allCNVs$size)
maxCNVSize <- max(allCNVs$size)

#Summary of deletions

totalDels <- dim(allDels)[1]
mean_nDels <- mean(as.numeric(sampleDetails$nDel))
std_devDels <- sd(as.numeric(sampleDetails$nDel))

Dels_below_two_std_dev <-sampleDetails[sampleDetails$nDel < mean_nDels-(2*std_devDels), ]
nDels_below_two_std_dev <- dim(Dels_below_two_std_dev)[1]

Dels_above_two_std_dev <-sampleDetails[sampleDetails$nDel > mean_nDels+(2*std_devDels), ]
nDels_above_two_std_dev <- dim(Dels_above_two_std_dev)[1]

meanDelSize <- mean(allDels$size)
minDelSize <-min(allDels$size)
maxDelSize <- max(allDels$size)

#Summary of Duplications
totalDups <- dim(allDups)[1]
mean_nDups <- mean(as.numeric(sampleDetails$nDup))
std_devDups <- sd(as.numeric(sampleDetails$nDup))

Dups_below_two_std_dev <-sampleDetails[sampleDetails$nDup <= mean_nDups-(2*std_devDups), ]
nDups_below_two_std_dev <- dim(Dups_below_two_std_dev)[1]

Dups_above_two_std_dev <-sampleDetails[sampleDetails$nDup >= mean_nDups+(2*std_devDups), ]
nDups_above_two_std_dev <- dim(Dups_above_two_std_dev)[1]


meanDupsize <- mean(allDups$size)
minDupSize <-min(allDups$size)
maxDupSize <- max(allDups$size)

allCNVranges <- GRanges(seqnames = allCNVs$chr,
                        ranges = IRanges(start = allCNVs$start,
                                         end = allCNVs$end),
                        type= allCNVs$type,
                        size= allCNVs$size,
                        sampleName= allCNVs$sampleName)
CNVolaps <- countOverlaps(query = allCNVranges,subject = allCNVranges)

#Stats Based on Cnv Size are below
"
#sampleDetails[i,] <-rbind(c(sampleName[i,1],nCNVs,
   #                           mean(CNV$size),
    #                         min(CNV$size),
     #                         max(CNV$size),
      #                       min(CNV[CNV$type=="<DEL>",6]),max(CNV[CNV$type=="<DEL>",6]),
       #                       nDups,mean(CNV[CNV$type=="<DUP>",6]),sd(CNV[CNV$type=="<DUP>",6]),
        ##                     mean(CNV$nTargetExons)))
  
#Summary of allCNVs
totalCNV <- dim(allCNVs)[1]
mean_nCNVs <-  mean(as.numeric(sampleDetails$nCNV)) 
std_devCNV <- sd(sampleDetails$nCNV)
meanCNVsize <- mean(allCNVs$size)
minCNVSize <- min(allCNVs$size)
maxCNVSize <- max(allCNVs$size)
CNV_below_two_std_dev <- allCNVs[allCNVs$size <= meanCNVsize-(2*std_devCNV), ]
CNV_above_two_std_dev <- allCNVs[allCNVs$size >= meanCNVsize+(2*std_devCNV), ]
nDels_above_two_std_dev <- dim(CNV_above_two_std_dev[CNV_above_two_std_dev$type=="<DEL>",])[1]
nDups_above_two_std_dev <- dim(CNV_above_two_std_dev[CNV_above_two_std_dev$type=="<DUP>",])[1]
nDels_below_two_std_dev <- dim(CNV_below_two_std_dev[CNV_below_two_std_dev$type=="<DEL>",])[1]
nDups_below_two_std_dev <- dim(CNV_below_two_std_dev[CNV_below_two_std_dev$type=="<DUP>",])[1]
#Summary of deletions
totalDels <- dim(allDels)[1]
meanDelSize <- mean(allDels$size)
std_devDels <- sd(allDels$size)
minDelSize <-min(allDels$size)
maxDelSize <- max(allDels$size)
Dels_below_two_std_dev <-allDels[allDels$size <= meanDelSize-(2*std_devDels), ]
nDels_below_two_std_dev <- dim(Dels_below_two_std_dev)[1]
Dels_above_two_std_dev <-allDels[allDels$size <= meanDelSize+(2*std_devDels), ]
nDels_above_two_std_dev <- dim(Dels_above_two_std_dev)[1]
#Summary of Duplications
totalDups <- dim(allDups)[1]
meanDupsize <- mean(allDups$size)
std_devDups <- sd(allDups$size)
minDupSize <-min(allDups$size)
maxDupSize <- max(allDups$size)
Dups_below_two_std_dev <-allDups[allDups$size <= meanDupsize-(2*std_devDups), ]
nDups_below_two_std_dev <- dim(Dups_below_two_std_dev)[1]
Dups_above_two_std_dev <-allDups[allDups$size <= meanDupsize+(2*std_devDups), ]
nDups_above_two_std_dev <- dim(Dups_above_two_std_dev)[1]
"




