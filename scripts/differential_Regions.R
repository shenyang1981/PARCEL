#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("data.table"));
suppressWarnings(suppressPackageStartupMessages(library(adagio)));

kadane = function(x,pos=NULL,minscore=0,dscore=0,poslist=NULL){
  if(!is.null(pos) & length(x) == length(pos)){
    tmp = rep(dscore,max(pos)-min(pos)+1);
    tmp[pos-min(pos)+1]=x;
    pos = min(pos):max(pos);
  } else {
    tmp = x;
    pos=1:length(x);
  }
  tmpsub = maxsub(tmp,inds=T,compiled=T);
  if(tmpsub$sum < minscore){
    if(is.null(poslist)){
      poslist[[1]] = list("start"=1,"end"=2,"maxScore"=-10);
    }
    return(poslist);
  } else {
    # if there is a sub-array with score passed cutoff, split whole array into two sub-arrays
    leftpos = tmpsub$inds[1]+min(pos)-1;
    rightpos = tmpsub$inds[2]+min(pos)-1;
    resultlist = list("start"=leftpos,"end"=rightpos,"maxScore"=tmpsub$sum);
    
    # store current max window
    if(is.null(poslist)){
      poslist=list();
      poslist[[1]] = resultlist;
    } else {
      poslist[[length(poslist)+1]]=resultlist;
    }
    
    # handle left subarray
    if(tmpsub$inds[1]>1){
      xleft = tmp[1:(tmpsub$inds[1]-1)];
      posleft = pos[1:(tmpsub$inds[1]-1)]; # real position of left subarray
      if(length(xleft)==1 & sum(xleft)>minscore){
        poslist[[length(poslist)+1]] = list("start"=posleft[1],"end"=posleft[1],"maxScore"=xleft);
      } else {
        poslist = kadane(xleft,pos=posleft,minscore=minscore,dscore=dscore,poslist=poslist);    
      }
    }
    
    # handle right subarray
    if(tmpsub$inds[2]<length(tmp)){
      xright = tmp[(tmpsub$inds[2]+1):length(tmp)];
      posright = pos[(tmpsub$inds[2]+1):length(tmp)]; # real position of left subarray
      if(length(xright)==1 & sum(xright)>minscore){
        poslist[[length(poslist)+1]] = list("start"=posright[1],"end"=posright[1],"maxScore"=xright);
      } else {
        poslist = kadane(xright,pos=posright,minscore=minscore,dscore=dscore,poslist=poslist);    
      }
    }
    
  }
  return(poslist);
}

evalue = function(n = NULL,sumscore=NULL){
  lambda = 0.862871
  K = 0.0809635
  X <- sumscore-log(n,base=exp(1))/lambda;
  Prob <- K*exp(-lambda*X);
  return(Prob);
}

lambda <- 0.862871
K <- 0.0809635

args = commandArgs(T);
args = args[-1];
print(args);
covcutoff = as.numeric(args[1]); # 10
conditions = args[2]; # treatment
resultdir = args[3]; 
covinfo = args[4]; # "covinfo_{treatment}.Rdata", 
load(covinfo);

print(conditions)
load(paste(resultdir,"etTable_",conditions,".Rdata",sep=""))
etTable = etTable[,score:=log(0.1)-log(PValue)];
#totalSite = dim(etTable)[1];
totalSite = sum(rowSums(v1all[,c(-1,-2),with=F]) > (NCOL(v1all)-2));

output = etTable[,do.call(rbind.data.frame,kadane(x=score,pos=pos,minscore=5,dscore=-10)),by=chr]
output = output[,c("E_value","pos","winSize","geneID","pstart","pend"):=
              list(evalue(n=totalSite,sumscore=maxScore),
                    round((start+end)/2,0),
                    end-start+1,
                    gsub("(.*):(.*)","\\1",chr),
                    start,
                    end)];
output = output[maxScore>=5,];
save(output,file=paste(resultdir,"fastq2_",conditions,"_output10.Rdata",sep=""))
print(paste("no. of output10",NROW(output),sep=":"));
