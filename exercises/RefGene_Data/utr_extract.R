codingTx <-read.csv(file ="C:/Users/Narayanan/Desktop/summer_worklog/coding_transcripts.csv",header = T,sep = ",")
fiveUTR <- data.frame()
threeUTR <- data.frame()
for (i in 1:dim(codingTx)[1]){
  
  exStart = strsplit(x = codingTx[i,10], split = ",")
  exEnd = strsplit(x= codingTx[i,11], split = "," )
  
  if (codingTx[i,4] == "+") {
    
     if (codingTx[i,5] == codingTx[i,7] && codingTx[i,7] != exStart[[1]][1]){ 
      
      fiveUTR[i,1:7] <- rbind(cbind("fiveUTR ",codingTx[i,5], exStart[[1]][1], codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,7] == exStart[[1]][1]){
      
      fiveUTR[i,1:7] <- rbind(cbind("fiveUTR ", codingTx[i,5], codingTx[i,7], codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,5] == exStart[[1]][1]){
      
      fiveUTR[i,1:7] <- rbind(cbind("fiveUTR ", codingTx[i,5], codingTx[i,7], codingTx[i,2:4], codingTx[i,13]))
    } 
    
    
     if (codingTx[i,6] == codingTx[i,8] && codingTx[i,8] != exEnd[[1]][codingTx[i,9]]){
      
      threeUTR[i,1:7] <- rbind(cbind("threeUTR ", exEnd[[1]][codingTx[i,9]], codingTx[i,5], codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if(codingTx[i,6] != codingTx[i,8] && codingTx[i,8] == exEnd[[1]][codingTx[i,9]]){
      
      threeUTR[i,1:7] <- rbind(cbind("threeUTR ",codingTx[i,8], codingTx[i,6], codingTx[i,2:4], codingTx[i,13]))
    }
    
    else if(codingTx[i,6] != codingTx[i,8] && codingTx[i,6] == exEnd[[1]][codingTx[i,9]]){
      
      threeUTR[i,1:7] <- rbind(cbind("threeUTR ",codingTx[i,8], codingTx[i,6], codingTx[i,2:4], codingTx[i,13]))
    }
  }
  
  if (codingTx[i,4] == "-"){
    
      if (codingTx[i,5] == codingTx[i,7] && codingTx[i,7] != exStart[[1]][1]){
      
      threeUTR[i,1:7] <- rbind(cbind("threeUTR ", codingTx[i,5],exStart[[1]][1], codingTx[i,2:4], codingTx[i,13]))
    }
    
    
    else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,7] == exStart[[1]][1]){
      
      threeUTR[i,1:7] <- rbind(cbind("threeUTR ",codingTx[i,5], codingTx[i,7], codingTx[i,2:4], codingTx[i,13]))
    }
    
    
    else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,5] == exStart[[1]][1]){
      
      threeUTR[i,1:7] <- rbind(cbind("threeUTR ",codingTx[i,5], codingTx[i,7], codingTx[i,2:4], codingTx[i,13]))
    }
    
     if (codingTx[i,6] == codingTx[i,8] && codingTx[i,8] != exEnd[[1]][codingTx[i,9]]){
      
      fiveUTR[i,1:7] <- rbind(cbind("fiveUTR ", exEnd[[1]][codingTx[i,9]], codingTx[i,5], codingTx[i,2:4], codingTx[i,13]))
    }
    else if(codingTx[i,6] != codingTx[i,8] && codingTx[i,8] == exEnd[[1]][codingTx[i,9]]){
      
      fiveUTR[i,1:7] <- rbind(cbind("fiveUTR ", codingTx[i,8], codingTx[i,6], codingTx[i,2:4], codingTx[i,13]))
    }
    else if(codingTx[i,6] != codingTx[i,8] && codingTx[i,6] == exEnd[[1]][codingTx[i,9]]){
      
      fiveUTR[i,1:7] <- rbind(cbind("fiveUTR ", codingTx[i,8], codingTx[i,6], codingTx[i,2:4], codingTx[i,13]))
    }
  }
}


colnames(fiveUTR) <- c("type","start","end","transcript_id","chrom","strand","geneName")
colnames(threeUTR) <- c("type","start","end","transcript_id","chrom","strand","geneName")

fiveUTR <- fiveUTR[complete.cases(fiveUTR$type),]
threeUTR <- threeUTR[complete.cases(threeUTR$type),]

write.csv(x=fiveUTR,file = "fiveUTR.csv",row.names = F)
write.csv(x=threeUTR,file = "threeUTR.csv",row.names = F)

all_utr <- merge(x = fiveUTR,y = threeUTR, by="chrom")
colnames(all_utr) <- c("chrom","strand","geneName","transcript_id","utr.1","start.1","end.1","utr.2","start.2","end.2")
write.csv(x = all_utr, file = "all_utr.csv", row.names = F)

