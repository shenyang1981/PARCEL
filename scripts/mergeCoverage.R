#!/usr/bin/env Rscript

# convert tab coverage into data.table: v1all
# row is chromsome:position and column is library id.

suppressPackageStartupMessages(library("data.table"))
args = commandArgs(T);
args = args[-1];

outdir = args[1]; 
dir.create(outdir);

covfile = args[2]; # allcov.wide.min20.txt.gz
samplefile = args[3]; #sampleList_Fly.txt
usedBatch = args[4]; #batch1
countinfo = args[5]; # combined_v1all.Rdata
ismerge = ifelse(!is.null(args[6]) & args[6]=="T",T,F);

sampleInfo = read.table(samplefile,header=T,sep="\t",stringsAsFactors = F);
sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == usedBatch,];

if(ismerge){
  sampleInfo[,"LibID"] = paste(sampleInfo[,"Condition"],sampleInfo[,"Replicates"],sep="__");
}
sampleInfo = unique(sampleInfo[,c("Condition","LibID")]);

v1all <- NULL

# load coverages
if(grepl("gz$",covfile)){
  iname = paste("zcat ",covfile,sep="");
}
v1all = fread(iname,header = T);
v1all = v1all[,c("V1","V2",sampleInfo[,"LibID"]),with=F];
setnames(v1all,1:2,c("chr","pos"));

# for transcriptome, we exclude reads mapped to negative strand of transcript
v1all = v1all[substr(chr,nchar(chr)-3,nchar(chr))!=":Neg",];
# exclude position zero
v1all = v1all[pos!=0,];
#save(v1all,file=paste(outdir,countinfo,sep=""));
save(v1all,file=countinfo);
