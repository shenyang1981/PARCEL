#!/usr/bin/env Rscript
source("/home/sheny/gseq/prog/parcel/scripts/sunScripts/useful_Rscript.R")
source("~/gseq/prog/parcel/rlibrary.R");
suppressPackageStartupMessages(library("data.table"));
#setwd("/mnt/projects/sunm/others/yeast/analysis2")

lambda <- 0.862871
K <- 0.0809635

args = commandArgs(T);
args = args[-1];
print(args);
covcutoff = as.numeric(args[1]); # 10
conditions = args[2]; # atp
resultdir = args[3]; # analysis_combined/
covinfo = args[4]; # "covinfo_met.Rdata", 
load(covinfo);

print(conditions)
load(paste(resultdir,"etTable_",conditions,".Rdata",sep=""))
etTable = etTable[,score:=log(0.1)-log(PValue)];
#totalSite = dim(etTable)[1];
totalSite = sum(rowSums(v1all[,c(-1,-2),with=F]) > (NCOL(v1all)-2));
output = etTable[,do.call(rbind.data.frame,kadaneShen(x=score,pos=pos,minscore=5,dscore=-10)),by=chr]
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