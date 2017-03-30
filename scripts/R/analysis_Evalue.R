library("reshape2");
library("data.table")
library("LSD");
library('RColorBrewer');
library("lattice");
library("statmod");
library("bedr");
library(polyester);
library("Hmisc");
library(rrcov);

#library(RSkittleBrewer);
#trop = RSkittleBrewer('tropical');
trop = c("darkorange","dodgerblue","hotpink","limegreen","yellow");
#library(genefilter)
#library(RUVSeq)

library(sva)

source("~/gseq/prog/rlib/common_lib.R");
source("~/gseq/prog/parcel/rlibrary.R");

setwd("~/gseq/prog/parcel/B.subtilis/");

#genomesize = 4215606;
# transcriptome size 3775323
# CDS size 3697728

# sample information





# coverage vs. sequence depth ---------------------------------------------



# coverage from RNA-seq

covinfo = read.table("PRJEB12568_TermSeq/coverages/allcov.txt.gz",header=F,stringsAsFactors = F);

colnames(covinfo) = c("chr","position","cov","sample");

covinfoW = dcast(covinfo,chr+position~sample,value.var = "cov");

covinfoW[,c(-1,-2)] = apply(covinfoW[,c(-1,-2)],2,function(x){x[is.na(x)]=0;return(x)});


covinfoM = read.table("PRJEB12568_TermSeq/coverages/allcov.merged.txt.gz",header=F,stringsAsFactors = F);

myprop = 10^(-1*(seq(0,3,0.1)))
covcutoff = c(5,10,20,40);
downS = data.frame("totalC"=0,"covcut"=0,"percentG"=0);
k = 1;
for(i in 1:length(myprop)){
  for(j in 1:length(covcutoff)){
    tmpx = downSample(x = covinfoM[,3],prop = myprop[i]);
    y = percentGenomeCov(x = tmpx,genomesize = genomesize,covcutoff=covcutoff[j]);
    downS[k,"totalC"] = sum(tmpx);
    downS[k,"covcut"] = covcutoff[j];  
    downS[k,"percentG"] = y;
    k = k+1;
  }
}

pdf("coverage_vs_sequence_depth_RNASeq.pdf");
mycols = c("black","red","blue","brown");
plot(0,0,xlim=range(downS[,"totalC"]),ylim=range(downS[,"percentG"]),pch='.',col="white",
     xlab="Total Mapped Reads",
     ylab="Proportion of Genome is Covered");
#plot(1e4,1e4,xlim=c(2e5,4e5),ylim=c(0,0.001),pch='.',col="white");
for(i in 1:length(covcutoff)){
  x = downS[downS[,"covcut"]==covcutoff[i],"totalC"];
  y = downS[downS[,"covcut"]==covcutoff[i],"percentG"];
  lines(x,y,type='l',lwd=1.5,col=mycols[i]);
  tmp = apply(covinfoW[,c(-1,-2)],2,function(x){percentGenomeCov(x = x,genomesize = genomesize , covcutoff=covcutoff[i])})
  x1 = colSums(covinfoW[,c(-1,-2)]);
  y1 = as.numeric(tmp);
  points(colSums(covinfoW[,c(-1,-2)]),as.numeric(tmp),pch=20,col=mycols[i]);
}
legend("topleft",legend=paste("coverage: ",covcutoff,sep=""),col = mycols,lty=rep(1,length(covcutoff)));
dev.off();


# coverage from RNA-seq

covinfo = read.table("PRJEB12568_TermSeq/coverages/allcov.txt.gz",header=F,stringsAsFactors = F);

colnames(covinfo) = c("chr","position","cov","sample");

covinfoW = dcast(covinfo,chr+position~sample,value.var = "cov");

covinfoW[,c(-1,-2)] = apply(covinfoW[,c(-1,-2)],2,function(x){x[is.na(x)]=0;return(x)});


covinfoM = read.table("PRJEB12568_TermSeq/coverages/allcov.merged.txt.gz",header=F,stringsAsFactors = F);



#binsize = 100;
mybins = c(50,100,150,200);
for(binsize in mybins){
  myprop = 10^(-1*(seq(0,3,0.1)))
  covcutoff = c(5,10,20,40);
  downS = data.frame("totalC"=0,"covcut"=0,"percentG"=0,"maxPercentG"=0);
  k = 1;
  for(i in 1:length(myprop)){
    for(j in 1:length(covcutoff)){
      tmpx = downSample(x = covinfoM[,3],prop = myprop[i]);
      y = percentGenomeCov(x = tmpx,genomesize = genomesize,covcutoff=covcutoff[j],bin=binsize,coordinates=covinfoM[,2]);
      y1 = percentGenomeCov(x = tmpx,genomesize = genomesize,covcutoff=covcutoff[j],bin=binsize,coordinates=covinfoM[,2],maxcov = T);
      downS[k,"totalC"] = sum(tmpx);
      downS[k,"covcut"] = covcutoff[j];  
      downS[k,"percentG"] = y;
      downS[k,"maxPercentG"] = y1; # coverage for maxcov in bin
      k = k+1;
    }
  }
  
  pdf(paste("coverage_vs_sequence_depth_bin_RNASeq",binsize,".pdf",sep=""));
  
  mycols = c("black","red","blue","brown");
  plot(0,0,xlim=range(downS[,"totalC"]),ylim=range(downS[,"percentG"]*1.5),pch='.',col="white",
       xlab="Total Mapped Reads",
       ylab="Proportion of Genome is Covered",main=paste("window size: ",binsize));
  #plot(1e4,1e4,xlim=c(2e5,4e5),ylim=c(0,0.001),pch='.',col="white");
  for(i in 1:length(covcutoff)){
    x = downS[downS[,"covcut"]==covcutoff[i],"totalC"];
    y = downS[downS[,"covcut"]==covcutoff[i],"percentG"];
    lines(x,y,type='l',lwd=1.5,col=mycols[i]);
    tmp = apply(covinfoW[,c(-1,-2)],2,function(x){percentGenomeCov(x = x,genomesize = genomesize , covcutoff=covcutoff[i],bin=binsize,coordinates = covinfoW[,2])})
    x1 = colSums(covinfoW[,c(-1,-2)]);
    y1 = as.numeric(tmp);
    points(colSums(covinfoW[,c(-1,-2)]),as.numeric(tmp),pch=20,col=mycols[i]);
  }
  legend("topleft",legend=paste("coverage: ",covcutoff,sep=""),col = mycols,lty=rep(1,length(covcutoff)),cex = 0.7);
  
  plot(0,0,xlim=range(downS[,"totalC"]),ylim=range(downS[,"maxPercentG"]*1.5),pch='.',col="white",
       xlab="Total Mapped Reads",
       ylab="Proportion of windows with at least one site passed coverage cutoff",main=paste("window size: ",binsize));
  #plot(1e4,1e4,xlim=c(2e5,4e5),ylim=c(0,0.001),pch='.',col="white");
  for(i in 1:length(covcutoff)){
    x = downS[downS[,"covcut"]==covcutoff[i],"totalC"];
    y = downS[downS[,"covcut"]==covcutoff[i],"maxPercentG"];
    lines(x,y,type='l',lwd=1.5,col=mycols[i]);
    tmp = apply(covinfoW[,c(-1,-2)],2,function(x){percentGenomeCov(x = x,genomesize = genomesize , covcutoff=covcutoff[i],bin=binsize,coordinates = covinfoW[,2],maxcov=T)})
    x1 = colSums(covinfoW[,c(-1,-2)]);
    y1 = as.numeric(tmp);
    points(colSums(covinfoW[,c(-1,-2)]),as.numeric(tmp),pch=20,col=mycols[i]);
  }
  legend("topleft",legend=paste("coverage: ",covcutoff,sep=""),col = mycols,lty=rep(1,length(covcutoff)),cex = 0.7);
  dev.off();
}


# quality control ---------------------------------------------------------

# adapter trimming and mapping statistics

trimInfo = read.table("trimAndMap/trimSummary.txt",header=F,sep="\t");
colnames(trimInfo) = c("LibID","Sequenced","AfterTrimmed","WithAdapter");

mapInfo = read.table("readsToGenome//mapSummary.txt",header=F,sep="\t");
colnames(mapInfo) = c("LibID","Input","Mapped","UniqMapped","Unmapped");

covinfoW = read.table("coverages/allcov.wide.min20.txt.gz",header=T,sep="\t",stringsAsFactors = F);

# PCA analysis

sampleInfo = read.table("sampleList_Bsub.txt",header=T,sep="\t",stringsAsFactors = F);

tmpdata = covinfoW[,sampleInfo[,"LibID"]];
colnames(tmpdata) = paste(sampleInfo[,"Condition"],sampleInfo[,"LibID"],sep="__");
tmpdata = tmpdata[rowSums(tmpdata)>100,]
tmpdata = apply(tmpdata,2,function(x){x=x/sum(x)*1e6;return(x)})
t_tmpdata = log2(t(tmpdata)+1);
pca_proc <- prcomp(t_tmpdata[,apply(t_tmpdata, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE,retX=TRUE)

# library(rgl)
# plot3d(pca_proc$x[,1:3],size=5)
# text3d(pca_proc$x[,1:3],texts=rownames(t_tmpdata),cex=0.7)

library("pca3d");
pca3d(pca_proc,group = sampleInfo[,"SeqBatch"],show.labels = rownames(t_tmpdata),radius = 1.5,labels.col = "grey");
snapshotPCA3d("pca3D_all.png");

library(rrcov)
pdf("robustPCA.pdf");
seqbatch = unique(sampleInfo[,"ExpermentalBatch"]);
outlierResult = NULL;
for(i in 1:length(seqbatch)){
  if(i ==1 ){
    tmpmatrix = t_tmpdata;
    pcaHub <- PcaHubert(tmpmatrix[,apply(tmpmatrix, 2, var, na.rm=TRUE) != 0])
    plot(pcaHub,pch=20,cex=1.5);
    outliers <- which(pcaHub@flag=='FALSE');
    if(length(outliers)!=0){
      tmp = data.frame("SeqBatch"="All","Outliers" = names(outliers));
      if(is.null(outlierResult)){
        outlierResult = tmp;
      } else {
        outlierResult = rbind(outlierResult,tmp);
      }
    } 
  }
  tmpmatrix = t_tmpdata[sampleInfo[,"ExpermentalBatch"]==seqbatch[i],]
  pcaHub <- PcaHubert(tmpmatrix[,apply(tmpmatrix, 2, var, na.rm=TRUE) != 0])
  plot(pcaHub,pch=20,cex=1.5);
  outliers <- which(pcaHub@flag=='FALSE');
  if(length(outliers)!=0){
    tmp = data.frame("ExpermentalBatch"=seqbatch[i],"Outliers" = names(outliers));
    if(is.null(outlierResult)){
      outlierResult = tmp;
    } else {
      outlierResult = rbind(outlierResult,tmp);
    }
  }
}
dev.off();
write.table(outlierResult,file="outliers.txt",col.names=T,row.names=F,sep="\t",quote=F);

allResult = merge(sampleInfo[,c("LibID","Condition","Replicates","SeqBatch","ExpermentalBatch")],trimInfo[,c("LibID","Sequenced","WithAdapter","AfterTrimmed")],by="LibID",all.x=T);

allResult = merge(allResult,mapInfo[,c("LibID","Mapped","UniqMapped")],by="LibID",all.x=T);

allResult[,"Outlier"] = "Good";

allResult[!is.na(match(allResult[,"LibID"],gsub("(.*)__(.*)","\\2",outlierResult[outlierResult[,"SeqBatch"]!="All",2]))),"Outlier"] = "Outlier";

allResultM = melt(allResult,id.vars = c("LibID","Condition","Outlier"),measure.vars = c("Sequenced","WithAdapter","AfterTrimmed","Mapped","UniqMapped"),variable.name = "Summary")

allResultM = allResultM[order(allResultM[,"Outlier"],allResultM[,""]),];
allResultM[,"value"] = log10(allResultM[,"value"]);
allResultM[,"ID"] = paste(allResultM[,"LibID"],allResultM[,"Condition"],sep="__");
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
pdf("SummaryOfProcessing.pdf",height=14,width=12);
barchart(ID~value|Summary,allResultM,groups =Outlier,
         xlab="Number of Reads (log10)",
         par.strip.text=list(cex=1.1,col="white"),
         par.settings=my.settings,box.ratio = 1, 
         panel = function(...) {   # Code for panel function by rcs (http://stackoverflow.com/questions/3220702/display-values-in-stacked-lattice-barchart-r) - thank you very much!
           panel.barchart(...);
           panel.grid(h=FALSE, v=-1);
           panel.abline(h=FALSE,v=log10(5e7),col="black",lty=2);
         },
         stack = F,layout = c(5,1),horizontal=T,
         auto.key = list(columns=2, space="top",                          
                         cex=0.8, size=1.4, adj=1,
                         between=0.2, between.colums=0.1, 
                         size = 1.3, points = FALSE, rectangles = TRUE))
dev.off();

write.table(allResult,file="processingSummary.xls",col.names=T,row.names=F,sep="\t",quote=F);


# Test differential Sites -------------------------------------------------


binsize = 1;

sampleInfo = read.table("sampleList_Bsub_buffered.txt",header=T,sep="\t",stringsAsFactors = F);

covinfoW = read.table("coverages/allcov.wide.min1.txt.gz",header=T,sep="\t",stringsAsFactors = F);
colnames(covinfoW)[1:2] = c("chr","pos");
# using SAM to test RiboSwitch

sampleInfo = read.table("sampleList_SAM.txt",header=T,sep="\t",stringsAsFactors = F);

cond1 = "control";
cond2 = "sam";

flagVsall = F;

tmpSampleInfo = sampleInfo[,c("Condition","LibID")];

if(flagVsall==T){
  tmpSampleInfo[tmpSampleInfo[,"Condition"]!=cond2,"Condition"] = cond1;  
}


mycountFull = covinfoW[,c("pos",tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"])];
mycount = mycountFull[rowSums(mycountFull[,-1])>5,];
mycountData = mycount[,-1];
rownames(mycountData) = mycount[,1];

coninfo = data.frame(
  "cond" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"Condition"],
  "cols" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"]);

nobatch = simulateData(counts = mycountData,isbatch = F,batchstr = 1,groupstr = 2,numGenes = 10000,propGrp = 0.2,propBatch = 0.8,myseed = 10000)

checkPCA(nobatch$data);

withbatch = simulateData(counts = mycountData,isbatch = T,batchstr = 1,groupstr = 2,numGenes = 10000,propGrp = 0.2,propBatch = 0.8,myseed = 10000,gcoeffs = nobatch$betamatrix[,1])

checkPCA(withbatch$data);


withbatchweak = simulateData(counts = mycountData,isbatch = T,batchstr = 0.2,groupstr = 2,numGenes = 10000,propGrp = 0.2,propBatch = 0.8,myseed = 10000,gcoeffs = nobatch$betamatrix[,1])

#tmp = cbind(nobatch$data,withbatch$data,withbatchweak$data);
colnames(nobatch$data) = c("Con_Rep1","Con_Rep2","Treat_Rep1","Treat_Rep2");
colnames(withbatch$data) = paste(c("Con_Rep1","Con_Rep2","Treat_Rep1","Treat_Rep2"),"_BAdded",sep="");
tmp = cbind(nobatch$data,withbatch$data);
pdf("Sim_Batch_Heatmap.pdf");
pheatmap(cor(nobatch$data,method='s'))
pheatmap(cor(withbatch$data,method='s'))
dev.off();
library("LSD");
heatpairs(log2(tmp+1))

# Test for no batch effects

mydeseq = diffSite(mycount = nobatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "DEseq",coninfo=coninfo);

mydeseq[,"g"] = nobatch$betamatrix[as.numeric(as.character(mydeseq[,"id"])),1];

myedgeR = diffSite(mycount = nobatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "edgeR",coninfo=coninfo,exactTest=T);

myedgeR[,"g"] = nobatch$betamatrix[as.numeric(as.character(myedgeR[,"id"])),1];

all = merge(mydeseq,myedgeR,by="id");
all[,"id"] = as.numeric(as.character(all[,"id"]))
table(all[,"id"]<2000,all[,"pvalue.x"]>log10(0.05)*-1);
table(all[,"id"]<2000,all[,"pvalue.y"]>log10(0.05)*-1)

cor(all[,"pvalue.x"],abs(all[,"g.x"]),method='s')
cor(all[,"pvalue.y"],abs(all[,"g.x"]),method='s')

# Test for batch effects

mydeseqB = diffSite(mycount = withbatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "DEseq",coninfo=coninfo);

mydeseqB[,"g"] = withbatch$betamatrix[as.numeric(as.character(mydeseqB[,"id"])),1];

myedgeRB = diffSite(mycount = withbatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "edgeR",coninfo=coninfo,exactTest=F);

myedgeRB[,"g"] = withbatch$betamatrix[as.numeric(as.character(myedgeRB[,"id"])),1];

allB = merge(mydeseqB,myedgeRB,by="id");
allB[,"id"] = as.numeric(as.character(all[,"id"]))
table(allB[,"id"]<2000,allB[,"pvalue.x"]>log10(0.05)*-1);
table(allB[,"id"]<2000,allB[,"pvalue.y"]>log10(0.05)*-1)

cor(allB[,"pvalue.x"],abs(allB[,"g.x"]),method='s')
cor(allB[,"pvalue.y"],abs(allB[,"g.x"]),method='s')

# Test for batch effects without batch effect model
coninfo[,"batch"] = factor(c(0,1,0,1));
mydeseqM = diffSite(mycount = nobatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "DEseq",coninfo=coninfo,batcheffect = T);

mydeseqM[,"g"] = nobatch$betamatrix[as.numeric(as.character(mydeseqM[,"id"])),1];

myedgeRM = diffSite(mycount = nobatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "edgeR",coninfo=coninfo,exactTest=F,batcheffect = T);

myedgeRM[,"g"] = nobatch$betamatrix[as.numeric(as.character(myedgeRM[,"id"])),1];

allBM = merge(mydeseqBM,myedgeRBM,by="id");
allBM[,"id"] = as.numeric(as.character(all[,"id"]))
table(allBM[,"id"]<2000,allBM[,"pvalue.x"]>log10(0.05)*-1);
table(allBM[,"id"]<2000,allBM[,"pvalue.y"]>log10(0.05)*-1)

cor(allBM[,"pvalue.x"],abs(allBM[,"g.x"]),method='s')
cor(allBM[,"pvalue.y"],abs(allBM[,"g.x"]),method='s')


# Test for batch effects with batch effect model
coninfo[,"batch"] = factor(c(0,1,0,1));
mydeseqBM = diffSite(mycount = withbatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "DEseq",coninfo=coninfo,batcheffect = T);

mydeseqBM[,"g"] = withbatch$betamatrix[as.numeric(as.character(mydeseqBM[,"id"])),1];

myedgeRBM = diffSite(mycount = withbatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "edgeR",coninfo=coninfo,exactTest=F,batcheffect = T);

myedgeRBM[,"g"] = withbatch$betamatrix[as.numeric(as.character(myedgeRBM[,"id"])),1];

allBM = merge(mydeseqBM,myedgeRBM,by="id");
allBM[,"id"] = as.numeric(as.character(all[,"id"]))
table(allBM[,"id"]<2000,allBM[,"pvalue.x"]>log10(0.05)*-1);
table(allBM[,"id"]<2000,allBM[,"pvalue.y"]>log10(0.05)*-1)

cor(allBM[,"pvalue.x"],abs(allBM[,"g.x"]),method='s')
cor(allBM[,"pvalue.y"],abs(allBM[,"g.x"]),method='s')

# comparison of Pvalue between nobatch vs. batch

comdiffResult = function(result1=NULL,result2=NULL,title=""){
  deseqcom = merge(result1,result2,by="id");
  deseqcom[,"id"] = as.numeric(as.character(deseqcom[,"id"]));
  myrange = quantile(c(deseqcom[,"pvalue.x"],deseqcom[,"pvalue.y"]),c(0.01,0.99))*1.3;
  plot(deseqcom[,"pvalue.x"],deseqcom[,"pvalue.y"],pch=20,cex=0.5,
       xlab="-log10(Pvalue) Batch Effect not been modelled",
       ylab="-log10(Pvalue) Batch Effect modelled",
       col=trop[2],xlim=myrange,ylim=myrange,main=title);
  points(deseqcom[deseqcom[,"id"]<2000,"pvalue.x"],deseqcom[deseqcom[,"id"]<2000,"pvalue.y"],col=trop[3],pch=20,cex=0.5);
  points(deseqcom[deseqcom[,"id"]>2000,"pvalue.x"],deseqcom[deseqcom[,"id"]>2000,"pvalue.y"],col=trop[2],pch=20,cex=0.5);
  abline(0,1,col=trop[1],lty=2);
  legend("bottomright",legend = c("Sig Diff Gene","Non-Sig Diff Gene"),col=c(trop[3],trop[2]),pch=c(20,20));
}
# deseq

#without batch effects

pdf("Sim_Pcom_NoBatch_DEseq.pdf");
comdiffResult(result1 =mydeseq,result2=mydeseqM,title="DESeq2");
dev.off();

pdf("Sim_Pcom_NoBatch_edgeR.pdf");
comdiffResult(result1 =myedgeR,result2=myedgeRM,title="EdgeR");
dev.off();

# with batch effects

pdf("Sim_Pcom_Batch_DEseq.pdf");
comdiffResult(result1 =mydeseqB,result2=mydeseqBM,title="DESeq2");
dev.off();

pdf("Sim_Pcom_Batch_edgeR.pdf");
comdiffResult(result1 =myedgeRB,result2=myedgeRBM,title="EdgeR");
dev.off();

# real SAM data

compareBatch = function(mycountData=NULL,coninfo=NULL,cond1=NULL,cond2=NULL,diffmethod="DEseq",prefix="",mincov=10,bmincov=40){
  
  realdeseq = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = mincov,method = diffmethod,coninfo=coninfo,autobatch = T,bmincov=bmincov);
  
  realdeseqNoB = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = mincov,method = diffmethod,coninfo=coninfo,autobatch = F,bmincov=bmincov);
  
  all = merge(realdeseq,realdeseqNoB,by="id");
  all[,"id"] = as.numeric(as.character(all[,"id"]));
  all[,"start"] = all[,"id"];
  all[,"end"] = all[,"id"]+1;
  all[,"chr"] = "chr";
  all[,"locus"] = bedToLocus(all[,c("chr","start","end")]);
  
  posR = read.table("sumResult/positiveSAM.bed",header=F,stringsAsFactors = F);
  posR[,2] = posR[,2]-100;
  posR[,3] = posR[,3]+100;
  
  png(paste(prefix,"_Manhattan.png",sep=""),width = 1200,height = 800,pointsize=20,res=72);
  plot(all[,"id"],all[,"pvalue.x"],
       ylim=range(c(all[,"pvalue.x"],all[,"pvalue.y"])),
       cex=0.2,col="grey",xlab="genomic coordinate",ylab="-log10(P)");
  abline(h=3,lty=2,lwd=0.3,col="grey");
  abline(v=posR[!grepl("Not",posR[,4]),2],col="grey",lwd=0.5)
  abline(v=posR[grepl("Not",posR[,4]),2],col=trop[4],lwd=0.5,lty=2)
  points(all[all[,"pvalue.x"]>3,"id"],all[all[,"pvalue.x"]>3,"pvalue.x"],cex=0.3,col=trop[1])
  points(all[all[,"pvalue.y"]>3,"id"],all[all[,"pvalue.y"]>3,"pvalue.y"],cex=0.3,col=trop[2])
  legend("topright",legend=c("Considering Batch","Not Considering Batch"),pch=c(20,20),col=c(trop[1],trop[2]));
  dev.off();
  
  png(paste(prefix,"_Manhattan_Batch.png",sep=""),width = 1200,height = 800,pointsize=20,res=72);
  plot(all[,"id"],all[,"pvalue.x"],
       ylim=range(c(all[,"pvalue.x"],all[,"pvalue.y"])),
       cex=0.2,col="grey",xlab="genomic coordinate",ylab="-log10(P)");
  abline(h=3,lty=2,lwd=0.3,col="grey");
  abline(v=posR[!grepl("Not",posR[,4]),2],col="grey",lwd=0.5)
  abline(v=posR[grepl("Not",posR[,4]),2],col=trop[4],lwd=0.5,lty=2)
  points(all[all[,"pvalue.x"]>3,"id"],all[all[,"pvalue.x"]>3,"pvalue.x"],cex=0.3,col=trop[1])
  dev.off();
  
  png(paste(prefix,"_Manhattan_NotBatch.png",sep=""),width = 1200,height = 800,pointsize=20,res=72);
  plot(all[,"id"],all[,"pvalue.x"],
       ylim=range(c(all[,"pvalue.x"],all[,"pvalue.y"])),
       cex=0.2,col="grey",xlab="genomic coordinate",ylab="-log10(P)");
  abline(h=3,lty=2,lwd=0.3,col="grey");
  abline(v=posR[!grepl("Not",posR[,4]),2],col="grey",lwd=0.5)
  abline(v=posR[grepl("Not",posR[,4]),2],col=trop[4],lwd=0.5,lty=2)
  points(all[all[,"pvalue.y"]>3,"id"],all[all[,"pvalue.y"]>3,"pvalue.y"],cex=0.3,col=trop[2])
  dev.off();
    
  posR.sort = bedr.sort.region(posR[!grepl("Not",posR[,4]),c(1:4)],check.chr=F);
  result.sort = bedr.sort.region(all[,c("chr","start","end")],check.chr=F);
  is.region <- in.region(result.sort, posR.sort,check.chr=F);
  overlapSite = result.sort[is.region,]
  
  allPos = all[!is.na(match(all[,"locus"],bedToLocus(overlapSite))),]
  plot(allPos[,3],allPos[,6],xlab="-log10(P) Considering Batch",ylab="-log10(P) Not Considering Batch",pch=20,col=trop[1],xlim=range(c(allPos[,3],allPos[,6])),ylim=range(c(allPos[,3],allPos[,6])),main="Sites within Known RiboSwitch Region",cex=0.7);
  abline(0,1,col=trop[4],lty=2);
  
  allNonePos = all[is.na(match(all[,"locus"],bedToLocus(overlapSite))),]
  smoothScatter(allNonePos[,3],allNonePos[,6],xlab="-log10(P) Considering Batch",ylab="-log10(P) Not Considering Batch",pch=20,col=trop[2],xlim=range(c(allNonePos[,3],allNonePos[,6])),ylim=range(c(allNonePos[,3],allNonePos[,6])),main="Sites Not in Known RiboSwitch Region",cex=0.5);
  abline(0,1,col=trop[4],lty=2);
  return(all);
}

# only use two libraries

cond1 = "control";
cond2 = "sam";

flagVsall = F;

tmpSampleInfo = sampleInfo[,c("Condition","LibID")];

if(flagVsall==T){
  tmpSampleInfo[tmpSampleInfo[,"Condition"]!=cond2,"Condition"] = cond1;  
}


mycountFull = covinfoW[,c("pos",tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"])];
mycount = mycountFull[rowSums(mycountFull[,-1])>5,];
mycountData = mycount[,-1];
rownames(mycountData) = mycount[,1];

coninfo = data.frame(
  "cond" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"Condition"],
  "cols" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"]);

pdf("Sim_Real_DESeq_Two.pdf");
compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="DEseq",prefix="Sim_Real_DESeq_Two",bmincov=0);
dev.off();

pdf("Sim_Real_edgeR_Two.pdf");
compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="edgeR",prefix="Sim_Real_edgeR_Two",bmincov=0);
dev.off();

flagVsall = T;

tmpSampleInfo = sampleInfo[,c("Condition","LibID")];

if(flagVsall==T){
  tmpSampleInfo[tmpSampleInfo[,"Condition"]!=cond2,"Condition"] = cond1;  
}

mycountFull = covinfoW[,c("pos",tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"])];
mycount = mycountFull[rowSums(mycountFull[,-1])>5,];
mycountData = mycount[,-1];
rownames(mycountData) = mycount[,1];

coninfo = data.frame(
  "cond" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"Condition"],
  "cols" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"]);

pdf("Sim_Real_DESeq_Three.pdf");
tmp = compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="DEseq",prefix="Sim_Real_DESeq_Three",bmincov=0);
dev.off();

pdf("Sim_Real_edgeR_Three.pdf");
tmp = compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="edgeR",prefix="Sim_Real_edgeR_Three",bmincov=0);
dev.off();

pdf("Sim_Real_DESeq_Three_40.pdf");
compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="DEseq",prefix="Sim_Real_DESeq_Three_40",bmincov=40);
dev.off();

pdf("Sim_Real_edgeR_Three_40.pdf");
compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="edgeR",prefix="Sim_Real_edgeR_Three_40",bmincov=40);
dev.off();

# poor all together

sampleInfoAll = read.table("sampleList_Bsub.txt",header=T,sep="\t",stringsAsFactors = F);


cond1 = "control";
cond2 = "sam";

flagVsall = T;

tmpSampleInfo = sampleInfoAll[,c("Condition","LibID")];

if(flagVsall==T){
  tmpSampleInfo[tmpSampleInfo[,"Condition"]!=cond2,"Condition"] = cond1;  
}


mycountFull = covinfoW[,c("pos",tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"])];
mycount = mycountFull[rowSums(mycountFull[,-1])>5,];
mycountData = mycount[,-1];
rownames(mycountData) = mycount[,1];

coninfo = data.frame(
  "cond" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"Condition"],
  "cols" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"],stringsAsFactor=F);

pdf("Sim_Real_DESeq_All.pdf");
compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="DEseq",prefix="Sim_Real_DESeq_All",bmincov=10,mincov=60);
dev.off();

pdf("Sim_Real_edgeR_All.pdf");
compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="edgeR",prefix="Sim_Real_edgeR_All",bmincov=10,mincov=60);
dev.off();


# test for removing those libraries with large batch effect
mincovtotal = 100;
tmpmatrix = t(mycountData[rowSums(mycountData)>mincovtotal,]);
pcaHub <- PcaHubert(tmpmatrix[,apply(tmpmatrix, 2, var, na.rm=TRUE) != 0],scale=T)
plot(pcaHub,pch=20,cex=1.5);
outliers <- which(pcaHub@flag=='FALSE');
if(length(outliers)!=0){
  coninfo = coninfo[is.na(match(coninfo[,"cols"],names(outliers))),];
}

# mincov = 50;
# tmp = mycountData;
# tmp[tmp<mincov]=NA;
# allpwcor = rcorr(as.matrix(tmp),type="spearman")$r;
pdf("Sim_Real_DESeq_Filter.pdf");
tmp = compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="DEseq",prefix="Sim_Real_DESeq_Filter",bmincov=10,mincov=60);
dev.off();

pdf("Sim_Real_edgeR_Filter.pdf");
tmp = compareBatch(mycountData=mycountData,coninfo=coninfo,cond1=cond1,cond2=cond2,diffmethod="edgeR",prefix="Sim_Real_edgeR_Filter",bmincov=10,mincov=60);
dev.off();



#diffmethod = "DEseq";
diffmethod = "edgeR";

realdeseq = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = diffmethod,coninfo=coninfo,autobatch = T);

realdeseqNoB = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = diffmethod,coninfo=coninfo,autobatch = F);

all = merge(realdeseq,realdeseqNoB,by="id");
all[,"id"] = as.numeric(as.character(all[,"id"]));
all[,"start"] = all[,"id"];
all[,"end"] = all[,"id"]+1;
all[,"chr"] = "chr";
all[,"locus"] = bedToLocus(all[,c("chr","start","end")]);

posR = read.table("sumResult/positiveSAM.bed",header=F,stringsAsFactors = F);
posR[,2] = posR[,2]-100;
posR[,3] = posR[,3]+100;


plot(all[,"id"],all[,"pvalue.x"],cex=0.2,col=trop[1])
points(all[,"id"],all[,"pvalue.y"],cex=0.2,col=trop[2])
abline(v=posR[!grepl("Not",posR[,4]),2],col="grey",lwd=0.5)
abline(v=posR[grepl("Not",posR[,4]),2],col=trop[4],lwd=0.5,lty=2)

posR.sort = bedr.sort.region(posR[,c(1:4)],check.chr=F);
result.sort = bedr.sort.region(all[,c("chr","start","end")],check.chr=F);
is.region <- in.region(result.sort, posR.sort,check.chr=F);
overlapSite = result.sort[is.region,]

allPos = all[!is.na(match(all[,"locus"],bedToLocus(overlapSite))),]
plot(allPos[,3],allPos[,6])
abline(0,1)

allNonePos = all[is.na(match(all[,"locus"],bedToLocus(overlapSite))),]
plot(allNonePos[,3],allNonePos[,6])
abline(0,1)

# Test for autobatch
mydeseqBAM = diffSite(mycount = withbatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "DEseq",coninfo=coninfo,batcheffect = T,autobatch = T);
mydeseqBAM[,"g"] = withbatch$betamatrix[as.numeric(as.character(mydeseqBAM[,"id"])),1];

allBMA = merge(mydeseqBM,mydeseqBAM,by="id");
allBMA[,"id"] = as.numeric(as.character(allBMA[,"id"]))
table(allBMA[,"id"]<2000,allBMA[,"pvalue.x"]>log10(0.05)*-1);
table(allBMA[,"id"]<2000,allBMA[,"pvalue.y"]>log10(0.05)*-1)

cor(allBMA[,"pvalue.x"],abs(allBMA[,"g.x"]),method='s')
cor(allBMA[,"pvalue.y"],abs(allBMA[,"g.x"]),method='s')

myedgeRBAM = diffSite(mycount = withbatch$data,cond1 = cond1,cond2 = cond2,maxrsum = 10,method = "edgeR",coninfo=coninfo,autobatch = T);
myedgeRBAM[,"g"] = withbatch$betamatrix[as.numeric(as.character(myedgeRBAM[,"id"])),1];

allBMA = merge(myedgeRBM,myedgeRBAM,by="id");
allBMA[,"id"] = as.numeric(as.character(allBMA[,"id"]))
table(allBMA[,"id"]<2000,allBMA[,"pvalue.x"]>log10(0.05)*-1);
table(allBMA[,"id"]<2000,allBMA[,"pvalue.y"]>log10(0.05)*-1)

cor(allBMA[,"pvalue.x"],abs(allBMA[,"g.x"]),method='s')
cor(allBMA[,"pvalue.y"],abs(allBMA[,"g.x"]),method='s')




# Simulate Batch effects --------------------------------------------------



# add batch effect

noiselevel = 0.5;
noisebase = 1.5;
#noise1=as.numeric(unlist(lapply(rowSums(mycountData),function(x){rbinom(1,x,noiselevel)})));
#noise2=as.numeric(unlist(lapply(rowSums(mycountData),function(x){rbinom(1,x,noiselevel)})));
noise1 = rbinom(dim(mycountData)[1],rowSums(mycountData),noiselevel);
noise2 = rbinom(dim(mycountData)[1],rowSums(mycountData),noiselevel);

noiseMatrix = covinfoW[rowSums(mycountFull[,-1])>5,c("RYY148","RYY149")];

noise1 = noiseMatrix[,1];
noise2 = noiseMatrix[,2];
mycountDataSim = mycountData;
mycountDataSim[,1]=mycountData[,1]+floor(rbinom(dim(mycountDataSim)[1],noise1,noiselevel)*noisebase);
mycountDataSim[,3]=mycountData[,3]+floor(rbinom(dim(mycountDataSim)[1],noise1,noiselevel)*noisebase);

checkPCA(mydata=mycountDataSim,mincov=10);
checkPCA(mydata=mycountData,mincov=10);

coninfo[,"batch"] = c("batch1","batch2","batch1","batch2");

mydeseq = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "DEseq",coninfo=coninfo);

myedgeR = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "edgeR",coninfo=coninfo);

myedgeRExact = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "edgeR",coninfo=coninfo,exactTest=T);

all = merge(myedgeR,myedgeRExact,by="id")

plot(all[,"pvalue.x"],all[,"pvalue.y"],col="red",pch=20,cex=0.5,xlim=c(0,40),ylim=c(0,40),main="GLM vs. Exact Test",xlab="GLM -log10(Pvalue)",ylab="Exact Test -log10(Pvalue)");
abline(0,1);

all = merge(myedgeR,mydeseq,by="id")

plot(all[,"pvalue.x"],all[,"pvalue.y"],col="red",pch=20,cex=0.5,xlim=c(0,40),ylim=c(0,40),main="edgeR GLM vs. DESeq GLM",xlab="edgeR GLM -log10(Pvalue)",ylab="DESeq GLM -log10(Pvalue)");
abline(0,1);

all = merge(myedgeRExact,mydeseq,by="id")

plot(all[,"pvalue.x"],all[,"pvalue.y"],col="red",pch=20,cex=0.5,xlim=c(0,40),ylim=c(0,40),main="edgeR Exact Test vs. DESeq GLM",xlab="edgeR Exact Test -log10(Pvalue)",ylab="DESeq GLM -log10(Pvalue)");
abline(0,1);


myedgeRSim = diffSite(mycount = mycountDataSim,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "edgeR",coninfo=coninfo);

myedgeRSimBatch = diffSite(mycount = mycountDataSim,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "edgeR",coninfo=coninfo,batcheffect=T);

mydeseqSim = diffSite(mycount = mycountDataSim,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "DEseq",coninfo=coninfo);

mydeseqSimBatch = diffSite(mycount = mycountDataSim,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "DEseq",coninfo=coninfo,batcheffect = T);

mydeseqBatch = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "DEseq",coninfo=coninfo,batcheffect = T);

all = merge(mydeseqBatch,mydeseqSimBatch,by="id")
plot(all[,"pvalue.x"],all[,"pvalue.y"],col="red",pch=20,cex=0.5,xlim=c(0,11),ylim=c(0,11),main="DEseq GLM vs. DESeq GLM (Batch)",xlab="DEseq GLM -log10(Pvalue)",ylab="DESeq GLM (Batch) -log10(Pvalue)");
abline(0,1);

all = merge(myedgeRSimBatch,mydeseqSimBatch,by="id")
plot(all[,"pvalue.x"],all[,"pvalue.y"],col="red",pch=20,cex=0.5,xlim=c(0,11),ylim=c(0,11),main="DEseq GLM vs. DESeq GLM (Batch)",xlab="DEseq GLM -log10(Pvalue)",ylab="DESeq GLM (Batch) -log10(Pvalue)");
abline(0,1);


#myedgeRSimBatch = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = 40,method = "edgeR",coninfo=coninfo,batcheffect=T);

all = merge(myedgeR,myedgeRSim,by="id")
all = merge(all,myedgeRSimBatch,by="id")
#plot(as.numeric(as.character(myedgeR[,"id"])),myedgeR[,"pvalue"])
plot(all[,"pvalue.x"],all[,"pvalue.y"],col="red",pch=20,cex=0.5,xlim=c(0,30),ylim=c(0,30));
abline(0,1);
points(all[,"pvalue.x"],all[,"pvalue"],col="blue",pch=20,cex=0.5);
abline(0,1);


# using Sun Miao's Data ---------------------------------------------------

load("SunMFastq/Bsub_6samples_v2all.Rdata") # v2all
v1all <- v2all

cond1 = "control";
cond2 = "sam";
mycount <- v1all[rowSums(v1all[,c(1,2,3,4,5,6)])>0,c(1,2,3,4,5,6)]
coninfo = data.frame(
  "cond" = c(rep("control",4),"sam","sam"),
  "cols" = colnames(mycount));

maxsum = 0;

mincov = 10;

sunedgeR = diffSite(mycount = mycount,cond1 = cond1,cond2 = cond2,maxrsum = maxsum,method = "edgeR",coninfo=coninfo,exactTest=T);

sunedgeR[,"id"] = as.numeric(as.character(sunedgeR[,"id"]));
sunedgeR[,"pos"] = as.numeric(as.character(sunedgeR[,"id"]));

#sunedgeR[sunedgeR[,"pos"]<4215606,"pos"] = sunedgeR[sunedgeR[,"pos"]<4215606,"pos"] + 2;
#sunedgeR[sunedgeR[,"pos"]>4215606,"pos"] = sunedgeR[sunedgeR[,"pos"]>4215606,"pos"] - 4215606 -2 ;
sunedgeR[,"score"] = log(0.1) - log(10^(-1*sunedgeR[,"pvalue"]));

# out raw pvalue to bed
forbed = aggregate(sunedgeR[sunedgeR[,"sumcount"]>mincov,"pvalue"],list(sunedgeR[sunedgeR[,"sumcount"]>mincov,"pos"]),max);
forbed[,"chr"] = "chr";
forbed[,"start"] = as.numeric(as.character(forbed[,1]));
forbed[,"end"] = forbed[,"start"]+1;
forbed[,"pvalue"] = 10^(-1*forbed[,2]);
write.table(forbed[,c("chr","start","end","pvalue")],file="sundataPvalue.bed",col.names=F,row.names=F,sep="\t",quote=F);

# get the region with maxiumal local score

maxScoreBin = do.call(rbind.data.frame,kadaneShen(x = sunedgeR[,"score"], pos = sunedgeR[,"pos"],minscore=5,dscore=-10));

# get the E value for each region
maxScoreBin = evalue(maxScoreBin = maxScoreBin,diffSiteResult = sunedgeR);


# Test for batch effect correction

sunedgeR = diffSite(mycount = mycount,cond1 = cond1,cond2 = cond2,maxrsum = maxsum,method = "edgeR",coninfo=coninfo,batcheffect = T);

sunedgeR[,"score"] = log(0.1) - log(10^(-1*sunedgeR[,"pvalue"]));
sunedgeR[,"id"] = as.numeric(as.character(sunedgeR[,"id"]));
sunedgeR[,"pos"] = as.numeric(as.character(sunedgeR[,"id"]));

# get the region with maxiumal local score

maxScoreBin = do.call(rbind.data.frame,kadaneShen(x = sunedgeR[,"score"], pos = sunedgeR[,"pos"],minscore=5,dscore=-10));

# get the E value for each region
maxScoreBin = evalue(maxScoreBin = maxScoreBin,diffSiteResult = sunedgeR);


maxScoreBin[,"chr"] = "chr";

maxScoreBin[,"start"] = ifelse(maxScoreBin[,"begin"]<4215606,maxScoreBin[,"begin"]+2,maxScoreBin[,"begin"]-4215606-2);
maxScoreBin[,"end"] = ifelse(maxScoreBin[,"end"]<4215606,maxScoreBin[,"end"]+2,maxScoreBin[,"end"]-4215606-2);

positeR = read.table("sumResult/positiveSAM.bed",header=F,sep="\t");
positeR = positeR[!grepl("Not",positeR[,4]),];
colnames(positeR)[1:3] = c("chr","start","end");
positeR[,"start"] = positeR[,"start"]-100;
positeR[,"end"] = positeR[,"end"]+100;


tmp = bedInterSect(region1=maxScoreBin,region2=positeR);
tmp = unique(tmp[,c(2:6)]);
load("SunMiaoAnalysis/Bsub2_v2all_SAM_output10.Rdata")


finalresult_miaodata_Deseq[finalresult_miaodata_Deseq[,"pos"]<4215606,"pos"] = finalresult_miaodata_Deseq[finalresult_miaodata_Deseq[,"pos"]<4215606,"pos"]+2;

finalresult_miaodata_Deseq[finalresult_miaodata_Deseq[,"pos"]>4215606,"pos"] = finalresult_miaodata_Deseq[finalresult_miaodata_Deseq[,"pos"]>4215606,"pos"] - 4215606 -2 ;

resultDeseq = aggregate(finalresult_miaodata_Deseq[,"pvalue"],list(finalresult_miaodata_Deseq[,"pos"]),max);
colnames(resultDeseq) = c("pos","pvalue");




allPvalue = merge(finalresultDeseq[,c("pos","pvalue")],resultEdgeR,by="pos")

allPvalue = merge(allPvalue,finalresult_mydata_EdgeR[,c("pos","pvalue")],by="pos");

allPvalue = merge(allPvalue,resultDeseq[,c("pos","pvalue")],by="pos");

colnames(allPvalue) = c("pos","Deseq_MyData","EdgeR_MiaoData","EdgeR_Mydata","Deseq_MiaoData")

allPvalue = allPvalue[rowSums(allPvalue[,-1])>7,];
pdf("comparison_Deseq_EdgeR.pdf");
heatpairs(as.matrix(allPvalue[,-1]),method='spearman');
dev.off();


finalresultEdgeR[finalresultEdgeR[,"pos"]<3000000 & finalresultEdgeR[,"pvalue"]>10,][1:10,]


save(etTable,file="Bsub_SAM_v2all_etTable.Rdata")

