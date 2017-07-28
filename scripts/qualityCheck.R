#!/usr/bin/env Rscript
library("edgeR");
library("DESeq2");
library("reshape2");
library("data.table")
library("LSD");
library('RColorBrewer');
library("lattice");
library("pheatmap");
library("plyr");


args = commandArgs(T);

# for debug
# args = c("coverages/allcov.wide.refgene.min1000.txt.gz","sampleList_shen_yang_Human.txt", "T", "10", "merged");
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


# merge sample

mergecols = function(x=NULL,mydata = NULL){
  x = as.character(x);
  if(length(x)==1){
    return(as.numeric(mydata[,x]));
  } else {
    return(rowSums(mydata[,x]));  
  }
}

mergerows = function(x=NULL,mydata = NULL){
  x = as.character(x);
  if(length(x)==1){
    return(as.numeric(mydata[x,]));
  } else {
    return(colSums(mydata[x,]));  
  }
}




# heatmap

#pheatmap(cor(tmpdata,method='s'),annotation_col = sampleInfo[,c("ExperimentalBatch","SeqBatch")],filename = "pairwiseHeatmap.pdf",width=16,height = 12)
library(Hmisc);

library(rrcov)
allResult = merge(sampleInfo[,c("LibID","Condition","Replicates","SeqBatch","ExperimentalBatch","ID")],trimInfo[,c("LibID","Sequenced","WithAdapter","AfterTrimmed")],by="LibID",all.x=T);
allResult = merge(allResult,mapInfo[,c("LibID","Mapped","UniqMapped")],by="LibID",all.x=T);  

allResult[,"Mappability"] = round(allResult[,"Mapped"]/allResult[,"AfterTrimmed"],4)*100;

allResult = allResult[order(allResult[,"ExperimentalBatch"],allResult[,"ID"]),];

allResultM = melt(allResult,id.vars = c("ID","Condition"),measure.vars = c("Sequenced","WithAdapter","AfterTrimmed","Mapped","UniqMapped","Mappability"),variable.name = "Summary")

allResultM = allResultM[order(allResultM[,"ID"]),];
allResultM[allResultM[,"Summary"]!="Mappability","value"] = log10(allResultM[allResultM[,"Summary"]!="Mappability","value"]);

myColours <- brewer.pal(6,"Blues");
my.settings <- list(
  #superpose.polygon=list(col=myColours[2:5], border="transparent"),
  superpose.polygon=list(col=c("red","blue","black"),border="transparent"),
  strip.background=list(col=myColours[6]),
  strip.border=list(col="black"),
  layout.heights = list(strip = 1.5),
  par.main.text =2,
  par.sub.text = 2
)
# pdf(summeryfile,height=figSize/1.2,width=figSize*1.5);
# barchart(ID~value|Summary,allResultM,
#          xlab="Number of Reads (log10)",
#          scales = list(x = list(relation = "free")),
#          par.strip.text=list(cex=1.1,col="white"),
#          par.settings=my.settings,box.ratio = 1, 
#          panel = function(...) {   # Code for panel function by rcs (http://stackoverflow.com/questions/3220702/display-values-in-stacked-lattice-barchart-r) - thank you very much!
#            panel.barchart(...);
#            panel.grid(h=FALSE, v=-1);
#            panel.abline(h=FALSE,v=log10(5e7),col="black",lty=2);
#          },
#          stack = F,layout = c(6,1),horizontal=T,
#          auto.key = list(columns=2, space="top",                          
#               cex=0.8, size=1.4, adj=1,
#               between=0.2, between.colums=0.1, 
#               size = 1.3, points = FALSE, rectangles = TRUE))
# dev.off();

write.table(allResult,file=processfile,col.names=T,row.names=F,sep="\t",quote=F);

save.image(file.path(outdir,"all.Rdata"));
