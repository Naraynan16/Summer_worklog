codingTx<-read.csv(file ="C:/Users/Narayanan/Desktop/summer_worklog/coding_transcripts.csv",header = T)

inStart <- list()
inEnd <- list() 
refgene_introns <- data.frame()
k=1
for (i in 1:dim(codingTx)[1]){
  l <- 1
  if (codingTx[i, 9] > 1){
    
      exStart = strsplit(x = codingTx[i, 10], split = ",")
      exEnd = strsplit(x= codingTx[i, 11], split = "," )
    
 for (j in 1:codingTx[i, 9]-1){
       
 
   if(isTRUE(exEnd[[1]][j] != exStart[[1]][j+1])){
    
    refgene_introns[k,1:7]<- cbind(codingTx[i,2:4],paste0("intron_",l),exEnd[[1]][j],exStart[[1]][j+1],codingTx[i,13])
      l <- l+1
      k <- k+1
   }
}
}
}
colnames(refgene_introns) <- c("transcript_id","chrom","strand","intron_no","start","end","geneName")
write.csv(x = refgene_introns,file = "refgene_introns.csv",row.names = F)
