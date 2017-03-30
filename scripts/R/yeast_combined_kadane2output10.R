source("/mnt/projects/sunm/bacteria/useful_Rscript.R")
#setwd("/mnt/projects/sunm/others/yeast/analysis2")

lambda <- 0.862871
K <- 0.0809635

args = commandArgs(T);
args = args[-1];
print(args);
rowinfotab = args[1]; # "/mnt/projects/sunm/others/yeast/fastq_combined/apt_1.tab"
covcutoff = as.numeric(args[2]); # 10
conditions = args[3]; # atp
resultdir = args[4]; # analysis_combined/
covinfo = args[5]; # "yeast_combined_v1all.Rdata", 
sffile = args[6]; # "yeast_combined_v1all_edgeR_sf.Rdata"
load(covinfo);

v1 <- readtab(rowinfotab)
l <- sapply(v1[1:length(v1)],length)+200
p1 <- as.character(unlist(mapply(rep,names(l),l)))
p2 <- as.character(unlist(lapply(l,function(x) c(1:x))))

for (i in conditions){
  print(i)
  load(paste(resultdir,"combined_v1all_etTable_",i,".Rdata",sep=""))
  score1 <- log(0.1)-log(etTable[[3]])
  score <- rep((-10),NROW(v1all))
  score[rowSums(v1all)>covcutoff] <- score1
  pos1 <- 1
  posN <- length(score)
  output <- kadane(pos1=pos1,posN=posN)
  outtemp <- cbind(output,label=1)
  count <- NROW(output)
  while (output[NROW(output),3] > 5){
    left <- kadane(pos1=pos1,posN=(output[NROW(output),1]-1))
    right <- kadane(pos1=(output[NROW(output),2]+1),posN=posN)
    outtemp <- rbind(outtemp,cbind(left,label=0),cbind(right,label=0))
    outtemp <- outtemp[order(outtemp[,1]),]  
    maxBegin <- outtemp[outtemp[,4]==0,1][which.max(outtemp[outtemp[,4]==0,3])]
    maxInd <- which(outtemp[,1]==maxBegin)
    output <- rbind(output,outtemp[maxInd,1:3])
    pos1 <- ifelse(maxInd>1,outtemp[maxInd-1,2]+1,1)
    if (maxInd<NROW(outtemp)){
      posN <- outtemp[maxInd+1,1]-1
    } else posN <- length(score)
    # posN <- ifelse(maxInd>=NROW(outtemp),outtemp[maxInd+1,1]-1,length(score)) 
    outtemp[maxInd,4] <- 1
    count <- count+1
  }
  passed <- which(rowSums(v1all)>NCOL(v1all))
  n <- length(passed)
  if (NROW(output)>0){
     Prob <- numeric()
     for (j in 1:NROW(output)){
     	 index <- output[j,1]:output[j,2]
    	 index <- index[index%in%passed]
    	 X <- sum(score[index])-log(n)/lambda
    	 Prob[j] <- K*exp(-lambda*X)}
     output <- cbind(data.frame(output),pos=round((output[,2]+output[,1])/2),winSize=output[,2]-output[,1]+1,E_value=Prob)
     output <- output[output[,6]<10,]
     if (NROW(output)>0){
     	output <- data.frame(output,geneID=p1[output[,1]],pstart=p2[output[,1]],pend=p2[output[,2]])
  	save(output,file=paste(resultdir,"fastq2_",i,"_output10.Rdata",sep=""))
  	print(paste("no. of output10",NROW(output),sep=":"))}}
}
