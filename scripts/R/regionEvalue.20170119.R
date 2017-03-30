#library("reshape2");
#library("data.table")
#library("LSD");
#library('RColorBrewer');
#library("lattice");
#library("statmod");
#library("bedr");
#library(polyester);
#library("Hmisc");
#library(rrcov);
#suppressWarnings(suppressPackageStartupMessages(library(reshape2)));

#library(RSkittleBrewer);
#trop = RSkittleBrewer('tropical');

# Load arguments ----------------------------------------------------------

suppressPackageStartupMessages(library("argparse"));
suppressPackageStartupMessages(library("qqman"))
suppressPackageStartupMessages(library("bedr"))
source("~/gseq/prog/rlib/common_lib.R");
source("~/gseq/prog/parcel/rlibrary.R");
trop = c("darkorange","dodgerblue","hotpink","limegreen","yellow");
#combppath = "~/gseq/prog/parcel/scripts/combp.sh";
#browser();
parser <- ArgumentParser(description='Identify sites under structral change using PARCEL');

# specify our desired options
# by default ArgumentParser will add an help option
#parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
#parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")

parser$add_argument("-i", "--infile", type="character", nargs=1, help="site coverage as Input", metavar="allcov.wide.min1.txt.gz",required=T);

parser$add_argument("-s", "--samplefile", type="character", nargs=1, help="sample information", metavar="sampleList_Bsub_buffered.txt",required=T);

parser$add_argument("-o", "--outdir", type="character", nargs=1, help="prefix of directory for all output", metavar="parcelResult",required=T);

parser$add_argument("-t", "--treatment", type="character", nargs =1, help="name of treatment condition in samplelist file", metavar="sam",required=T);

parser$add_argument("-c", "--control", type="character", nargs =1, help="name of control condition in samplelist file", metavar="control",required=T);

parser$add_argument("-b", "--usedBatch", type="character", nargs =1, help="batch ID for comparison (defined in comparison batch column) samplelist file", metavar="control",required=T);

parser$add_argument("--mincov", default=20, type="double", metavar=20,help="sites with minimum coverage of will be used in test [default %(default)s]");

parser$add_argument("--downSampleProp", default=1, type="double", metavar=20,help="downsample 0.xx of original reads[default %(default)s]");

parser$add_argument("--nthread", default=4, type="double", metavar=4,help="number of threads for parallel computing (multiple chromsomes) [default %(default)s]");

parser$add_argument("--batchcov", default=20, type="double", metavar=20,help="sites with minimum coverage of will be used in estimating batch effect [default %(default)s]");

parser$add_argument("--seedlogP", default=3, type="double", metavar=3,help="only sites with seed pvalue will be used for estimating Evalue [default %(default)s]");

parser$add_argument("--seedCombp", default=3, type="double", metavar=3,help="only sites with seed pvalue will be used for combined-p [default %(default)s]");

parser$add_argument("--Ecutoff", default=1, type="double", metavar=1,help="Evalue cutoff for significant region [default %(default)s]");

parser$add_argument("--topCov", default=0.01, type="double", metavar=0.01,help="mask x friction of top coverage sites [default %(default)s]");

parser$add_argument("--ChrCov", default=0.001, type="double", metavar=0.001,help="drop chromosome if less than 0.xxx of chromosome is covered [default %(default)s]");

parser$add_argument("--mydescore", default=-10, type="double", metavar=0.01,help="mask x friction of top coverage sites [default %(default)s]");

parser$add_argument("--vsall", default="T", type="character", metavar="T",help="whether use other samples in sample list as control [default %(default)s]");

parser$add_argument("--posfile", default="", type="character", metavar="",help="bed file for positive control [default %(default)s]");

parser$add_argument("--chrsize", default="", type="character", metavar="",help="length of chromosome [default %(default)s]");

parser$add_argument("--batch", default="F", type="character", metavar="F",help="whether model batch effect [default %(default)s]");

parser$add_argument("--exactTest", default="F", type="character", metavar="F",help="whether use EdgeR exact Test [default %(default)s]");

parser$add_argument("--autobatch", default="T", type="character", metavar="T",help="whether detect batch effect automatically [default %(default)s]");

parser$add_argument("--autofilter", default="F", type="character", metavar="F",help="whether filter outlier using rosubst PCA [default %(default)s]");

parser$add_argument("--runcombp", default="F", type="character", metavar="F",help="whether run combined-p [default %(default)s]");

parser$add_argument("--splitByChr", default="F", type="character", metavar="F",help="whether run Test for each chromosome independently [default %(default)s]");

parser$add_argument("--method", default="edgeR", type="character", metavar="edgeR",help="package for testing differential sites [default %(default)s]");

parser$add_argument("--diffGene", default="F", type="character", metavar="F",help="whether test differential expressed Gene [default %(default)s]");

args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
#if ( args$verbose ) {
#  write("writing some verbose output to standard error...\n", stderr())
#}


covfile = args$infile;
outdir = args$outdir;
samplefile = args$samplefile;
usedBatch = args$usedBatch;
posfile = args$posfile;
chrsize = args$chrsize;
cond2 = args$treatment;
cond1 = args$control;

mincov = args$mincov;
minbatchcov = args$batchcov;
seedlogP = args$seedlogP;
seedCombp = args$seedCombp;
Ecutoff = args$Ecutoff;
topCov = 1-args$topCov;
chrCov = args$ChrCov;
method = args$method;
flagVsall=as.logical(args$vsall);
isbatch = as.logical(args$batch);
isexactTest = as.logical(args$exactTest);
autobatch = as.logical(args$autobatch);
autofilter = as.logical(args$autofilter);
runcombp = as.logical(args$runcombp);
splitByChr = as.logical(args$splitByChr);
diffGene = as.logical(args$diffGene);
nthread = args$nthread;
mydescore = args$mydescore;
downSampleProp = args$downSampleProp;

# load positive control file ----------------------------------------------
posRegion=NULL;
if(posfile!=""){
  posRegion = read.table(posfile,header=F,sep="\t",quote='"',stringsAsFactor=F);
  posRegion = bedr.sort.region(posRegion,check.chr=F,verbose = F);
}

# args = commandArgs(T);
# args = c("sampleList_Bsub_buffered.txt","coverages/allcov.wide.min1.txt.gz","control","sam","T","T","edgeR");
# samplefile = args[1];
# covfile = args[2];
# cond1 = args[3];
# cond2 = args[4];
# mincov = args[5];
# flagVsall=as.logical(args[6]);
# isbatch = as.logical(args[7]);
# method=args[8];



#Rscript ../regionEvalue.R "-i coverages/allcov.wide.min20.txt.gz -s sampleList_SAM_Three.txt -o samResult -t sam -c control"

# covfile = "coverages/allcov.wide.min20.txt.gz";
# outdir = "samResult";
# samplefile = "sampleList_SAM_Three.txt";
# 
# cond2 = "sam";
# cond1 = "control";
# 
# mincov = 20;
# minbatchcov = 20;
# seedlogP = 5;
# Ecutoff = 1;
# method = "edgeR";
# flagVsall=T;
# isbatch = T;
# isexactTest = F;
# autobatch = T;
# autofilter = T;
# splitByChr = F;
# nthread = 1;


# Prepare Reads Count -----------------------------------------------------

sampleInfo = read.table(samplefile,header=T,sep="\t",stringsAsFactors = F);

sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == usedBatch,];

covinfoW = read.table(covfile,header=T,sep="\t",stringsAsFactors = F,quote='"');
colnames(covinfoW)[1:2] = c("chr","pos");
covinfoW[is.na(covinfoW)]=0;
tmpSampleInfo = sampleInfo[,c("Condition","LibID")];
#browser();
if(flagVsall==T){
  tmpSampleInfo[tmpSampleInfo[,"Condition"]!=cond2,"Condition"] = cond1;  
}

mycountFull = covinfoW[,c("chr","pos",tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"])];

# downsampling

if(downSampleProp!=1){
  #mycountFull[,c(-1,-2)] = apply(mycountFull,2,function(x){downSample(x=x,prop = downSampleProp);});
}

#browser();

if(diffGene==F){
  covPerChr = data.frame(table(gsub("(.*):.*","\\1",mycountFull[,"chr"])));
  colnames(covPerChr) = c("chr","sites");
  #print(covPerChr);
  
  if(file.exists(chrsize)){
    chrlength = read.table(chrsize,header=F,sep="\t");
  } else {
    chrlength = aggregate(covinfoW[,"pos"],list(gsub("(.*):.*","\\1",covinfoW[,"chr"])),max);
  }
  colnames(chrlength) = c("chr","totallen");
  
  covPerChr = merge(covPerChr,chrlength,by="chr");
  excludeChrs = covPerChr[ covPerChr[,"sites"]/covPerChr[,"totallen"] < chrCov,"chr"];
  totalcov = rowSums(mycountFull[,c(-1,-2)]);
  mycount = mycountFull[rowSums(mycountFull[,c(-1,-2)])>mincov & totalcov<=quantile(totalcov,topCov) & is.na(match(gsub("(.*):.*","\\1",mycountFull[,"chr"]),excludeChrs)),];
} else {
  totalcov = rowSums(mycountFull[,c(-1,-2)]);
  mycount = mycountFull[rowSums(mycountFull[,c(-1,-2)])>mincov & totalcov<=quantile(totalcov,topCov),];
}

#print(dim(mycountFull)[1]);
#print(dim(mycount)[1]);
print(paste("Number of chromosomes retained: ",length(unique(gsub("(.*):.*","\\1",mycount[,"chr"]))),sep=""));

# Creat Directory for Result -----------------------------------------------------------

batchtag = ifelse(isbatch | autobatch,"Batch","NoBatch");
vsalltag = ifelse(flagVsall,"vsAll","onlyControl");
#batchtag = ifelse(isbatch,"Batch","NoBatch");
outdir = paste(outdir,usedBatch,cond1,cond2,batchtag,vsalltag,mincov,sep="_");
dir.create(outdir, showWarnings = F, recursive = T, mode = "0755");
setwd(outdir);

# record parameters
sink("parameters.txt", append=FALSE, split=FALSE)
params = sort(names(args));
for(i in 1:length(params)){
  cat(paste(params[i],": ",args[[params[i]]],"\n",sep=""));
}
sink();

cat(paste("Number of Sites are retained: ",dim(mycount)[1],"\n",sep=""));


cat("Start Test...\n");
# Run Test ----------------------------------------------------------------
#browser();
runTest = function(mycount=NULL,mincov=10,cond1=NULL,cond2=NULL,method="edgeR",tmpSampleInfo=NULL,isexactTest=F,bmincov=10,isbatch=T,autobatch=T,autofilter=T,seedlogP=3,plottofile=NULL,diffGene=F){
  if(!is.null(plottofile)){
    pdf(plottofile);
  }  
  
  mycountData = mycount[,c(-1,-2)];
  rownames(mycountData) = paste(mycount[,1],mycount[,2],sep=":");
  
  coninfo = data.frame(
    "cond" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"Condition"],
    "cols" = tmpSampleInfo[!is.na(match(tmpSampleInfo[,"Condition"],c(cond1,cond2))),"LibID"]);
  
  
  if(autofilter==T & dim(mycountData)[2]>4){ # filter outliters for sample size > 4
    suppressPackageStartupMessages(library("rrcov"))
    tmpmatrix = t(mycountData[rowSums(mycountData)>bmincov,]);
    pcaHub <- PcaHubert(tmpmatrix[,apply(tmpmatrix, 2, var, na.rm=TRUE) != 0],scale=T)
    plot(pcaHub,pch=20,cex=1.5,main="Outliers Robust PCA");
    outliers <- which(pcaHub@flag=='FALSE');
    if(length(outliers)!=0){
      coninfo = coninfo[is.na(match(coninfo[,"cols"],names(outliers))),];
      sink("outliers.txt", append=FALSE, split=FALSE)
      cat(paste(names(outliers),collapse="\n"));
      sink();
    }
  }
  
  allresult = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = mincov,method = method,coninfo=coninfo,isexactTest=isexactTest,bmincov=minbatchcov,batcheffect = isbatch,autobatch = autobatch);
  #browser();
  if(diffGene==T){
    diffResult=allresult$result;
    normresult = data.frame("id" = diffResult[,"id"],allresult$norcpm);
    rawresult = data.frame("id" = diffResult[,"id"],allresult$rawcount);
    #colnames(normresult)[1] = "#chr";
    #colnames(rawresult)[1] = "#chr";
    myresult = list("result"=diffResult,"norcpm"=normresult,"rawcount"=rawresult);
    #return(allresult);
    return(myresult);
  }
  diffResult=allresult$result;
  diffResult[,"chr"] = gsub("(.*):(.*)","\\1",diffResult[,"id"]);
  diffResult[,"pos"] = as.numeric(gsub("(.*):(.*)","\\2",diffResult[,"id"]));
  diffResult[,"score"] = log(0.1) - log(10^(-1*diffResult[,"pvalue"]));
  
  diffResult[,"bin"] = "No";
  allchr = unique(diffResult[,"chr"]);
  totalBin = NULL;
  #browser();
  for(mychr in allchr){
    maxScoreBin = do.call(rbind.data.frame,kadaneShen(x = diffResult[diffResult[,"chr"]==mychr,"score"], pos = diffResult[diffResult[,"chr"]==mychr,"pos"],minscore=seedlogP,dscore=mydescore));  
    
    # attach bin information to diffResult
    for(i in 1:dim(maxScoreBin)[1]){
      tmpbin = as.numeric(maxScoreBin[i,]);
      diffResult[diffResult[,"chr"]==mychr & diffResult[,"pos"]>=tmpbin[1] & diffResult[,"pos"]<=tmpbin[2],"bin"] = paste(mychr,":",format(tmpbin[1],scientific=F),"-",format(tmpbin[2],scientific=F),sep="")
    }
    if(is.null(totalBin)){
      totalBin = maxScoreBin;
    } else {
      totalBin = rbind(totalBin,maxScoreBin);
    }
  }
  if(!is.null(plottofile)){
    dev.off();
  }
  siteid = rownames(allresult$norcpm);
  sitechr = gsub("(.*):(.*)","\\1",siteid);
  sitepos = as.numeric(gsub("(.*):(.*)","\\2",siteid));
  mychrs = gsub("(.*):(.*)","\\1",sitechr);
  mystrands = gsub("(.*):(.*)","\\2",sitechr);
  mystrands[mystrands=="Pos"]="+";
  mystrands[mystrands=="Neg"]="-";
  normresult = data.frame("chr"=mychrs,"start"=sitepos,"end"=sitepos+1,"name"="site","score"=diffResult[,"pvalue"],"strand"=mystrands,"minchange"=allresult$minchange,allresult$norcpm);
  rawresult = data.frame("chr"=mychrs,"start"=sitepos,"end"=sitepos+1,"name"="site","score"=diffResult[,"pvalue"],"strand"=mystrands,"minchange"=allresult$minchange,allresult$rawcount);
  colnames(normresult)[1] = "#chr";
  colnames(rawresult)[1] = "#chr";
  myresult = list("result"=diffResult,"norcpm"=normresult,"rawcount"=rawresult);
  return(myresult);
}



result = NULL;
if(splitByChr == T){
  allchrs = unique(covinfoW[,"chr"]);
  if(nthread>1 & length(allchrs)>1){
    suppressPackageStartupMessages(library("doMC"));
    #cl <- makeCluster(nthread, type="SOCK")
    #clusterExport(cl, c("interBioTypeMatrix"))
    registerDoMC(nthread);
    
    cat(paste("Parallel computing on ",nthread," cores ...\n",sep=""));
    
    
    # Using rbind is very time consumming
    #result = foreach(icount(trail),.combine = rbind)%dopar%{
    
    result = foreach(i=1:length(allchrs),.combine = rbind)%dopar%{
      tmpcov = mycount[mycount[,"chr"]==allchrs[i],];
      tmp = runTest(mycount=tmpcov,mincov=mincov,cond1=cond1,cond2=cond2,method=method,tmpSampleInfo=tmpSampleInfo,isexactTest=isexactTest,bmincov=minbatchcov,isbatch=isbatch,autobatch=autobatch,autofilter=autofilter,seedlogP=seedlogP,plottofile=paste(allchrs[i],".pdf",sep=""),diffGene=diffGene);
      tmp$rawcount[,c(2,3)] = format(tmp$rawcount[,c(2,3)],scientific = F);
      tmp$norcpm[,c(2,3)] = format(tmp$norcpm[,c(2,3)],scientific = F);
      write.table(tmp$norcpm,file=paste(allchrs[i],".cov.txt",sep=""),col.names=colnames(tmp$norcpm),row.names=F,sep="\t",quote=F);
      write.table(tmp$rawcount,file=paste(allchrs[i],".cov.raw.txt",sep=""),col.names=colnames(tmp$norcpm),row.names=F,sep="\t",quote=F);
      tmp$result;
        
    }
    cat("Parallel computing is end\n");
    #registerDoSEQ();
    #stopCluster(cl);
  } else {
    for(i in 1:length(allchrs)){
      tmpcov = mycount[mycount[,"chr"]==allchrs[i],];
      tmp = runTest(mycount=tmpcov,mincov=mincov,cond1=cond1,cond2=cond2,method=method,tmpSampleInfo=tmpSampleInfo,isexactTest=isexactTest,bmincov=minbatchcov,isbatch=isbatch,autobatch=autobatch,autofilter=autofilter,seedlogP=seedlogP,plottofile=paste(allchrs[i],".pdf",sep=""),diffGene=diffGene);
      if(is.null(result)){
        result = tmp$result;
      } else {
        result = rbind(result,tmp$result);
      }
      tmp$rawcount[,c(2,3)] = format(tmp$rawcount[,c(2,3)],scientific = F);
      tmp$norcpm[,c(2,3)] = format(tmp$norcpm[,c(2,3)],scientific = F);
      write.table(tmp$norcpm,file=paste(allchrs[i],".cov.txt",sep=""),col.names=colnames(tmp$norcpm),row.names=F,sep="\t",quote=F);
      write.table(tmp$rawcount,file=paste(allchrs[i],".cov.raw.txt",sep=""),col.names=colnames(tmp$norcpm),row.names=F,sep="\t",quote=F);  
    }
  }
} else {
  tmp = runTest(mycount=mycount,mincov=mincov,cond1=cond1,cond2=cond2,method=method,tmpSampleInfo=tmpSampleInfo,isexactTest=isexactTest,bmincov=minbatchcov,isbatch=isbatch,autobatch=autobatch,autofilter=autofilter,seedlogP=seedlogP,plottofile="genome.pdf",diffGene=diffGene);
  write.table(tmp$norcpm,file="all.cov.txt",col.names=colnames(tmp$norcpm),row.names=F,sep="\t",quote=F);  
  write.table(tmp$rawcount,file="all.cov.raw.txt",col.names=colnames(tmp$norcpm),row.names=F,sep="\t",quote=F);  
  result = tmp$result;
}

if(diffGene==T){
  write.table(result,file="geneDiff.txt",col.names=T,row.names=F,sep="\t",quote=F);
  cat("Finish testing differential genes\n");
  quit("no");
}

#browser();
totalSite = dim(result)[1];
tmp = aggregate(result[result[,"bin"]!="No","score"],list(result[result[,"bin"]!="No","bin"]),function(x){evalue(n=totalSite,sumscore=sum(x))});
colnames(tmp) = c("bin","Evalue");
tmp[is.infinite(tmp[,2]),"Evalue"]=1e2;
result = merge(tmp,result,by="bin");
result[,"realchr"] = gsub("(.*):(.*)","\\1",result[,"chr"]);
result[,"strand"] = gsub("(.*):(.*)","\\2",result[,"chr"]);
result[result[,"strand"]=="Neg","strand"]="-";
result[result[,"strand"]=="Pos","strand"]="+";
result[,"start"] = result[,"pos"];
result[,"end"] = result[,"pos"]+1;
result[,"P"] = 10^(-1*result[,"pvalue"]);

cat("Output results...\n");

# Manhattan Plot ----------------------------------------------------------

#manfile = paste(outdir,"/manhattan.png",sep="");
manfile = "manhattan.png";
tmp = result[,c("realchr","pos","P")];
colnames(tmp) = c("CHR","BP","P");
tmp[,"SNP"] = paste(tmp[,"CHR"],":",tmp[,"BP"],sep="");
sigSite = tmp[result[,"Evalue"]<Ecutoff,"SNP"];

tmp[,"CHR"] = factor(tmp[,"CHR"],levels = unique(tmp[,"CHR"]));
chrlabs = as.character(unique(tmp[,"CHR"]));
tmp[,"CHR"] = as.numeric(tmp[,"CHR"]);

png(manfile,width = 1200,height = 800,pointsize=20,res=72);
manhattan(tmp,chrlabs=chrlabs,highlight = sigSite,
          suggestiveline = F,
          cex = 0.5, cex.axis = 0.8,col = c("grey", trop[1]))
dev.off();
options(scipen=0);

# Out put raw result ------------------------------------------------------
#browser();
#diffResultFile = paste(outdir,"/diffSites.txt.gz",sep="");
diffResultFile = "diffSites.txt.gz";
gz1 <- gzfile(diffResultFile, "w");
write.table(result,gz1,col.names=T,row.names=F,sep="\t",quote=F);
close(gz1);

# bed file

#bedResultFile = paste(outdir,"/diffSites.Pvalue.bed",sep="");
bedResultFile = "diffSites.Pvalue.bed";
write.table(result[,c("chr","start","end","P")],file=bedResultFile,col.names=F,row.names=F,sep="\t",quote=F);


# Significant Region
#sigFile = paste(outdir,"/sigRegion.bed",sep="");
sigFile = "sigRegion.bed";
sigRegion = unique(result[result[,"Evalue"]<Ecutoff,c("bin","Evalue")]);
if(dim(sigRegion)[1]!=0){
  sigRegion[,"chr"] = gsub("(.*):(.*):(.*)-(.*)","\\1",sigRegion[,"bin"],perl=T);
  sigRegion[,"chr"] = gsub("-","___",sigRegion[,"chr"]);
  sigRegion[,"strand"] = gsub("(.*):(.*):(.*)-(.*)","\\2",sigRegion[,"bin"],perl=T);
  sigRegion[sigRegion[,"strand"]=="Pos","strand"]="+";
  sigRegion[sigRegion[,"strand"]=="Neg","strand"]="-";
  sigRegion[,"start"] = as.numeric(gsub("(.*):(.*):(.*)-(.*)","\\3",sigRegion[,"bin"],perl=T));
  sigRegion[,"end"] = as.numeric(gsub("(.*):(.*):(.*)-(.*)","\\4",sigRegion[,"bin"],perl=T));
  
  sigRegion.sort = bedr.sort.region(sigRegion[,c("chr","start","end","bin","Evalue","strand")],check.chr=F,verbose=F);
  sigRegion.sort.merge = bedr(
    engine = "bedtools", 
    input = list(i = sigRegion.sort), 
    method = "merge", 
    params = "-s -d 100 -c 5 -o min",
    check.zero.based=F,
    check.chr=F,
    verbose=F
  );
  colnames(sigRegion.sort.merge) = c("chr","start","end","strand","Evalue");
  sigRegion.sort.merge[,"chr"] = gsub("___","-",sigRegion.sort.merge[,"chr"]);
  sigRegion.sort.merge[,"bin"] = bedToLocus(sigRegion.sort.merge[,c("chr","start","end")])
} else {
  sigRegion.sort.merge = data.frame("chr"="chr","start"=1,"end"=1,"bin"="NoSigRegion","Evalue"=100,"strand"="+");
}
write.table(unique(sigRegion.sort.merge[,c("chr","start","end","bin","Evalue","strand")]),file=sigFile,col.names=F,row.names=F,sep="\t",quote=F);

# run comb-p
seedCombpValue = 10^(-1*seedCombp);
if(runcombp==T){
  combpCmd = paste("~/gseq/prog/parcel/scripts/combp.sh diffSites.Pvalue.bed diffSites.Pvalue.sort.bed ",seedCombpValue," diff_1e",seedCombp,sep="");
  system(combpCmd);
}



#

# merge combin-p and Evalue Result ----------------------------------------

# library("bedr")

minPbyMethod = function(x){
  method = gsub("(.*):(.*)","\\1",unlist(strsplit(x,split=',')));
  pvalue = as.numeric(gsub("(.*):(.*)","\\2",unlist(strsplit(x,split=','))));
  result = list("combineP"=NA,"Evalue"=NA);
  if(length(unique(method))==1){
    method = unique(method);
    result[[method]]=min(pvalue);
  } else {
    #browser();
    result = tapply(pvalue,method,min);  
  }
  result[['combineP']] = ifelse(sum('combineP'==names(result))==0,NA,result[['combineP']]);
  result[['Evalue']] = ifelse(sum('Evalue'==names(result))==0,NA,result[['Evalue']]);
  return(result);
}

#browser();
combpfile = paste("diff_1e",seedCombp,".regions-t.bed",sep="");
siginfo = file.info(sigFile);
siginfo$size = ifelse(is.na(siginfo$size),0,siginfo$size);
combpinfo = file.info(combpfile);
combpinfo$size = ifelse(is.na(combpinfo$size),0,combpinfo$size);
if(siginfo$size>30 & combpinfo$size>50){
  a = read.table(combpfile,header=F,sep="\t",stringsAsFactors = F)
  a[,1] = gsub("(.*):(.*)","\\1",a[,1]);
  a[,1] = gsub("-","___",a[,1]);
  b = read.table(sigFile,header=F,sep="\t",stringsAsFactors = F)
  b[,1] = gsub("-","___",b[,1]);
  a.sort = bedr.sort.region(a,check.chr=F,verbose=F)
  b.sort = bedr.sort.region(b,check.chr=F,verbose=F)
  a.sort[,"id"]=paste("combineP",a.sort[,7],sep=":")
  b.sort[,"id"] = paste("Evalue",b[,5],sep=":");
  all = bedr.sort.region(rbind(a.sort[,c("V1","V2","V3","id")],b.sort[,c("V1","V2","V3","id")]),check.chr=F,verbose=F);
  all.merge = bedr(
    engine = "bedtools", 
    input = list(i = all), 
    method = "merge", 
    params = "-d 50 -c 4 -o distinct",
    check.zero.based=F,
    check.chr=F,
    verbose=F
  );
  if(!is.null(posRegion)){
    tmp = bedr.join.region(all.merge[,1:4],posRegion,check.chr=F,verbose = F);
    all.merge = tmp[,c(1,2,3,4,8)];
    colnames(all.merge)[5]="Positive";
  }
  all.merge[,c("combineP","Evalue")] = as.data.frame(do.call(rbind,lapply(all.merge[,4],function(x){unlist(minPbyMethod(x=x))})));
  colnames(all.merge)[1:3] = c("chr","start","end");
  #print(class(all.merge));
  all.merge[,1]=gsub("___","-",all.merge[,1]);
  write.table(unique(all.merge),file="mergeSigRegion.txt",col.names=T,row.names=F,sep="\t",quote=F);
}

cat("Finish testing differential sites\n");
quit("no");


