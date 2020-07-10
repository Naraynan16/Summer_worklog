codingTx<-read.csv(file ="coding_transcripts.csv",header = T)
for (i in 1:dim(codingTx)[1])
    {
    exStart = strsplit(x = codingTx[i,10], split = ",")
  exEnd = strsplit(x= codingTx[i,11], split = "," )
  
if (codingTx[i,5] == codingTx[i,7] && codingTx[i,5] == exStart[[1]][1])
    { 
      if (codingTx[i,4] == "+") 
        {
          print(cbind("5UTR nil", codingTx[i,13]))
      }
      else 
      {
        print(cbind("3UTR nil", codingTx[i,13]))
    }
  }
else if (codingTx[i,5] != codingTx[i,7] && codingTx[i,5] == exStart[[1]][1])
{
  if (codingTx[i,4] == "+") 
  {
    print(cbind("5UTR ",codingTx[i,5], exStart[[1]][1], codingTx[i,13]))
  }
  else 
  {
    print(cbind("3UTR ", codingTx[i,5],exStart[[1]][1], codingTx[i,13]))
  }
}
  
else if(codingTx[i,5] != codingTx[i,7] && codingTx[i,5] != exStart[[1]][1])
{
  if (codingTx[i,4] == "+") 
  {
    print(cbind("5UTR ", codingTx[i,5],  codingTx[i,7], codingTx[i,13]))
  }
  else  
  {
    print(cbind("3UTR ",codingTx[i,5],  codingTx[i,7], codingTx[i,13]))
  }
}
}