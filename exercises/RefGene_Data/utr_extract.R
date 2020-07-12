codingTx<-read.csv(file ="C:/Users/Narayanan/Desktop/summer_worklog/coding_transcripts.csv",header = T)
utr1 <- data.frame()
utr2 <- data.frame()
for (i in 1:dim(codingTx)[1])
{
  exStart = strsplit(x = codingTx[i,10], split = ",")
  exEnd = strsplit(x= codingTx[i,11], split = "," )
  
  if (codingTx[i,4] == "+") {
    
    if (codingTx[i,5] == codingTx[i,7] && codingTx[i,7] == exStart[[1]][1]){         
      
      utr1[i, ] <- rbind(cbind("5UTR ", "-", "-", codingTx[i,2:4], codingTx[i,13]))
      
    }
    
    else if (codingTx[i,5] == codingTx[i,7] && codingTx[i,7] != exStart[[1]][1]){ 
      
      utr1[i,1:7] <- rbind(cbind("5UTR ",codingTx[i,5], exStart[[1]][1], codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,7] == exStart[[1]][1]){
      
      utr1[i,1:7] <- rbind(cbind("5UTR ", codingTx[i,5], codingTx[i,7], codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,5] == exStart[[1]][1]){
      
      utr1[i,1:7] <- rbind(cbind("5UTR ", codingTx[i,5], codingTx[i,7], codingTx[i,2:4], codingTx[i,13]))
    } 
    
    if (codingTx[i,6] == codingTx[i,8] && codingTx[i,8] == exEnd[[1]][codingTx[i,9]]){
      
      utr2[i, ] <- rbind(cbind("3UTR ", "-", "-", codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if (codingTx[i,6] == codingTx[i,8] && codingTx[i,8] != exEnd[[1]][codingTx[i,9]]){
      
      utr2[i,1:7] <- rbind(cbind("3UTR ", exEnd[[1]][codingTx[i,9]], codingTx[i,5], codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if(codingTx[i,6] != codingTx[i,8] && codingTx[i,8] == exEnd[[1]][codingTx[i,9]]){
      
      utr2[i,1:7] <- rbind(cbind("3UTR ",codingTx[i,8], codingTx[i,6], codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if(codingTx[i,6] != codingTx[i,8] && codingTx[i,6] == exEnd[[1]][codingTx[i,9]]){
      
      utr2[i,1:7] <- rbind(cbind("3UTR ",codingTx[i,8], codingTx[i,6], codingTx[i,2:4], codingTx[i,13]))
    }
  }
  
  if (codingTx[i,4] == "-"){
    
    if (codingTx[i,5] == codingTx[i,7] && codingTx[i,7] == exStart[[1]][1]){ 
      
      utr1[i, ] <- rbind(cbind("3UTR ", "-", "-", codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if (codingTx[i,5] == codingTx[i,7] && codingTx[i,7] != exStart[[1]][1]){
      
      utr1[i,1:7] <- rbind(cbind("3UTR ", codingTx[i,5],exStart[[1]][1], codingTx[i,2:4], codingTx[i,13]))
    }
    
    
    else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,7] == exStart[[1]][1]){
      
      utr1[i,1:7] <- rbind(cbind("3UTR ",codingTx[i,5], codingTx[i,7], codingTx[i,2:4], codingTx[i,13]))
    }
    
    
    else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,5] == exStart[[1]][1]){
      
      utr1[i,1:7] <- rbind(cbind("3UTR ",codingTx[i,5], codingTx[i,7], codingTx[i,2:4], codingTx[i,13]))
    }
    
    if (codingTx[i,6] == codingTx[i,8] && codingTx[i,8] == exEnd[[1]][codingTx[i,9]]){
      
      utr2[i, ] <- rbind(cbind("5UTR ", "-", "-", codingTx[i,2:4], codingTx[i,13]))
    }  
    else if (codingTx[i,6] == codingTx[i,8] && codingTx[i,8] != exEnd[[1]][codingTx[i,9]]){
      
      utr2[i,1:7] <- rbind(cbind("5UTR ", exEnd[[1]][codingTx[i,9]], codingTx[i,5], codingTx[i,2:4], codingTx[i,13]))
    }
    else if(codingTx[i,6] != codingTx[i,8] && codingTx[i,8] == exEnd[[1]][codingTx[i,9]]){
      
      utr2[i,1:7] <- rbind(cbind("5UTR ", codingTx[i,8], codingTx[i,6], codingTx[i,2:4], codingTx[i,13]))
    }
    else if(codingTx[i,6] != codingTx[i,8] && codingTx[i,6] == exEnd[[1]][codingTx[i,9]]){
      
      utr2[i,1:7] <- rbind(cbind("5UTR ", codingTx[i,8], codingTx[i,6], codingTx[i,2:4], codingTx[i,13]))
    }
  }
}

colnames(utr1) <- c("type","start","end","transcript_id","chrom","strand","geneName")
colnames(utr2) <- c("type","start","end","transcript_id","chrom","strand","geneName")
all_utr <- merge(x = utr1,y = utr2, by="chrom")
colnames(all_utr) <- c("chrom","strand","geneName","transcript_id","utr.1","start.1","end.1","utr.2","start.2","end.2")
write.csv(x = all_utr, file = "all_utr.csv", row.names = F)
