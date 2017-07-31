#!/usr/bin/env Rscript

args = commandArgs(T);

sampleinfofile = args[1];
trimsum = args[2];
mapsum = args[3];
combatch = args[4];
outdir = ifelse(is.null(args[5]),"./",args[5]);
dir.create(outdir, showWarnings = F, recursive = T, mode = "0755");

summeryfile = paste(outdir,"SummaryOfProcessing.pdf",sep="/");

processfile = paste(outdir,"processingSummary.xls",sep="/");

trimInfo = read.table(trimsum,header=F,sep="\t");
colnames(trimInfo) = c("LibID","Sequenced","AfterTrimmed","WithAdapter");

mapInfo = read.table(mapsum,header=F,sep="\t");
colnames(mapInfo) = c("LibID","Input","Mapped","UniqMapped","Unmapped");

sampleInfo = read.table(sampleinfofile,header=T,sep="\t",stringsAsFactors = F);
if(combatch!="All"){
  sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"]==combatch,];  
}

sampleInfo[,"ExperimentalBatch"] = factor(sampleInfo[,"ExperimentalBatch"]);

sampleInfo[,"ID"] = paste(sampleInfo[,"Condition"],sampleInfo[,"LibID"],sep="__");
rownames(sampleInfo) = paste(sampleInfo[,"Condition"],sampleInfo[,"LibID"],sep="__");

allResult = merge(sampleInfo[,c("LibID","Condition","Replicates","SeqBatch","ExperimentalBatch","ID")],trimInfo[,c("LibID","Sequenced","WithAdapter","AfterTrimmed")],by="LibID",all.x=T);
allResult = merge(allResult,mapInfo[,c("LibID","Mapped","UniqMapped")],by="LibID",all.x=T);  

allResult[,"Mappability"] = round(allResult[,"Mapped"]/allResult[,"AfterTrimmed"],4)*100;

allResult = allResult[order(allResult[,"ExperimentalBatch"],allResult[,"ID"]),];

write.table(allResult,file=processfile,col.names=T,row.names=F,sep="\t",quote=F);

save.image(file.path(outdir,"all.Rdata"));
