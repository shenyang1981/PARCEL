suppressPackageStartupMessages(library("argparse"));
suppressPackageStartupMessages(library("bedr"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("shape"))
suppressPackageStartupMessages(library("wordcloud"))
source("~/gseq/prog/rlib/common_lib.R");
source("~/gseq/prog/parcel/rlibrary.R");

trop = c("darkorange","dodgerblue","hotpink","limegreen","brown");
#combppath = "~/gseq/prog/parcel/scripts/combp.sh";
#browser();
parser <- ArgumentParser(description='Plot the coverage for significant region');

# specify our desired options
# by default ArgumentParser will add an help option
#parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
#parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")


parser$add_argument("-o", "--outdir", type="character", nargs=1, help="prefix of directory for all output", metavar="coveragePlot",required=T);

parser$add_argument("-i", "--infile", type="character", nargs=1, help="significant regions as Input", metavar="sigRegionFinal.xls",required=T);

parser$add_argument("-s", "--samplefile", type="character", nargs=1, help="sample information", metavar="sampleList_Bsub_buffered.txt",required=T);

parser$add_argument("-b", "--usedBatch", type="character", nargs =1, help="batch ID for comparison (defined in comparison batch column) samplelist file", metavar="control",required=T);

parser$add_argument("-t", "--treatment", type="character", nargs =1, help="name of treatment condition in samplelist file", metavar="sam",required=T);

parser$add_argument("-c", "--control", type="character", nargs =1, help="name of control condition in samplelist file", metavar="control",required=T);


parser$add_argument("--covfile", type="character", nargs=1, help="bed file with coverage information", metavar="all.cov.txt",required=T);
parser$add_argument("--vsall", default="T", type="character", metavar="T",help="whether use other samples in sample list as control [default %(default)s]");
parser$add_argument("--replot", default="F", type="character", metavar="F",help="only replot coverage using covinfo.xls [default %(default)s]");
parser$add_argument("--extWinSize", default=200, type="double", metavar=20,help="window size for calculating the wFC[default %(default)s]");
parser$add_argument("--genome", default="F", type="character", metavar="hg19.fas",help="genome Fasta [default %(default)s]");
parser$add_argument("--getFasta", default="F", type="character", metavar="F",help="whether extract sequence from Genome [default %(default)s]");
parser$add_argument("--CDSBed", default="", type="character", metavar="annotation.bed",help="CDS annotation in BED format[default %(default)s]");

# parser$add_argument("-s", "--sitefile", type="character", default="", help="fold change and pvalue for every site", metavar="diffSites.txt.gz");
# parser$add_argument("--Ecutoff", default=1, type="double", metavar=20,help="Evalue cutoff[default %(default)s]");
# parser$add_argument("--wFCcutoff", default=0.25, type="double", metavar=20,help="weighted Fold Change cutoff[default %(default)s]");
# parser$add_argument("--extWinSize", default=100, type="double", metavar=20,help="window size for calculating the wFC[default %(default)s]");
# parser$add_argument("--beta", default=2, type="double", metavar=20,help="beta in cutoff function[default %(default)s]");



args <- parser$parse_args()

infile = args$infile;
outdir = args$outdir;
covfile = args$covfile;
samplefile = args$samplefile;
usedBatch = args$usedBatch;
cond2 = args$treatment;
cond1 = args$control;
flagVsall=as.logical(args$vsall);
ext.win.size = args$extWinSize;
replot = as.logical(args$replot);
getFasta = as.logical(args$getFasta);
genome = args$genome;
CDSBed=args$CDSBed;
# infile="testevalueFilter_1_0.25/sigRegionFinal.xls"
# covfile="testevalue_batch1_control_sam_Batch_vsAll_20/all.cov.txt"
# samplefile = "sampleList_Bsub.txt";
# usedBatch = "batch1";
# cond2 = "sam";
# cond1 = "control";
# outdir = "testout";
# flagVsall=T;
# ext.win.size=100;

sampleInfo = read.table(samplefile,header=T,sep="\t",stringsAsFactors = F);

sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == usedBatch,];

tmpSampleInfo = sampleInfo[,c("Condition","LibID")];

if(flagVsall==T){
  tmpSampleInfo[tmpSampleInfo[,"Condition"]!=cond2,"Condition"] = cond1;  
}

treatmentCols = tmpSampleInfo[tmpSampleInfo[,"Condition"]==cond2,"LibID"];
controlCols = tmpSampleInfo[tmpSampleInfo[,"Condition"]==cond1,"LibID"];

# get CDS regions
cds.sort = NULL;
if(file.exists(CDSBed)){
  cds = read.table(CDSBed,header=F,sep="\t",quote='"',stringsAsFactor=F);
  cds.sort = bedr.sort.region(cds,check.chr=F,check.zero.based = F,verbose = F);
}

dir.create(outdir, showWarnings = F, recursive = T, mode = "0755");

if(replot==F){
  sigregion = read.table(infile,header=T,sep="\t",stringsAsFactors = F)
  sigregion.bed = sigregion[,c("chr","start","end","region","evalue","strand","gene","weightFold")]
  sigMidpoint = round((sigregion.bed[,2] + sigregion.bed[,3])/2,0);
  sigWidth = sigregion.bed[,3] - sigregion.bed[,2] + 1;
  sigregion.bed[sigWidth<ext.win.size*2,2] = sigMidpoint[sigWidth<ext.win.size*2] -ext.win.size;
  sigregion.bed[sigWidth<ext.win.size*2,3] = sigMidpoint[sigWidth<ext.win.size*2] +ext.win.size;
  
  
  
  sigregion.sort = bedr.sort.region(sigregion.bed,check.chr=F,check.zero.based = F,verbose = F)
  
  sites = read.table(covfile,header=T,sep="\t",quote='"',comment.char='$',stringsAsFactors=F)
  sites.sort = bedr.sort.region(sites,check.chr=F,check.zero.based = F,verbose = F,check.merge = F)
  
  #sigsites = bedr.join.region(sigregion.sort,sites.sort,check.chr=F,check.merge = F,check.zero.based = F,verbose=F)
  sigsites = bedr(
    engine = "bedtools", 
    input = list(a = sigregion.sort ,b=sites.sort), 
    method = "intersect", 
    params = "-s -loj -sorted",
    check.chr=F,
    verbose=F
  );
  if(!is.null(cds.sort)){
    tmpsig = sigsites[,c("X.chr.b","start.b","end.b","name","score","strand")];
    tmpsig[,2] = as.numeric(tmpsig[,2]);
    tmpsig[,3] = as.numeric(tmpsig[,3]);
    tmpsig.sort = bedr.sort.region(tmpsig,check.chr=F,check.zero.based = F,verbose = F,check.merge = F)
    tmpCDSsites = bedr(
      engine = "bedtools", 
      input = list(a = tmpsig.sort  ,b=cds.sort[,1:6]), 
      method = "intersect", 
      params = "-s -loj -sorted",
      check.chr=F,
      verbose=F
    );
    sigsitesIndex = bedToLocus(sigsites[,c("X.chr.b","start.b","end.b")]);
    cdssiteIndex = bedToLocus(tmpCDSsites[tmpCDSsites[,"V4"]!=".",c("X.chr.b","start.b","end.b")]);
    sigsites[,"isCDS"] = "No";
    sigsites[!is.na(match(sigsitesIndex,cdssiteIndex)),"isCDS"] = "Yes";
  }
  #browser();
  sigsites[,"evalue"] = as.numeric(sigsites[,"evalue"]);
  sigsites[,"weightFold"] = as.numeric(sigsites[,"weightFold"]);
  sigsites[,c(treatmentCols,controlCols)] = apply(sigsites[,c(treatmentCols,controlCols)],2,as.numeric);
  
  # shift to 5' position
  sigsites[,"start.b"] = as.numeric(sigsites[,"start.b"]);
  sigsites[,"end.b"] = as.numeric(sigsites[,"end.b"]);
  sigsites[sigsites[,"strand"]=="-",c("start.b","end.b")] = sigsites[sigsites[,"strand"]=="-",c("start.b","end.b")]+2;
  sigsites[sigsites[,"strand"]=="+",c("start.b","end.b")] = sigsites[sigsites[,"strand"]=="+",c("start.b","end.b")]-1;
  
  setwd(outdir);
} else {
  setwd(outdir);
  sigsites = read.table("covinfo.xls",header=T,sep="\t",quote='"',stringsAsFactor=F);
}

# set outdir
#outdir = paste(outdir,ecut,foldcut,sep="_");
#browser();
ddply(sigsites,.(region),plotCoverage,treatmentCols=treatmentCols,controlCols=controlCols,genome=genome,getFasta=getFasta,cds.sort=cds.sort);

write.table(sigsites[,c("region","evalue","strand","gene","weightFold","start.b","score",treatmentCols,controlCols)],file="covinfo.xls",col.names=T,row.names=F,sep="\t",quote=F);

if(!is.null(cds.sort)){
  sigMidpoint = round((sigregion.sort[,2] + sigregion.sort[,3])/2,0);
  sigregion.sort[,2] =  sigMidpoint;
  sigregion.sort[,3] =  sigMidpoint+1;
  #browser();
  regionType = NULL;
  #first see whether it located in a upstream region of any CDS within 500bp
  upregion = bedr(
      engine = "bedtools", 
      input = list(a = sigregion.sort[,1:6] ,b=cds.sort[,1:6]), 
      method = "closest", 
      params = "-s -D b -d -io -id",
      check.chr=F,
      verbose=F
    );
  upregion = upregion[upregion[,"V10"]!=".",];
  upregion[,"V13"] = as.numeric(upregion[,"V13"]);
  if(sum(abs(upregion[,"V13"])<500)>0){
    upregions = upregion[abs(upregion[,"V13"])<500,c("V4","V10","V13")];
    colnames(upregions) = c("region","regionType","closestCDSDistance");
    upregions[,"regionType"] = paste("5UTR__",upregions[,"regionType"],sep="");
    regionType = upregions;
  }
  
  #exclude region from 5'UTR
  sigregion.sort.new = sigregion.sort[is.na(match(sigregion.sort[,4],regionType[,1])),1:6]
  if(dim(sigregion.sort.new)[1]!=0){
    downregion = bedr(
      engine = "bedtools", 
      input = list(a = sigregion.sort.new ,b=cds.sort[,1:6]), 
      method = "closest", 
      params = "-s -D b -d -iu",
      check.chr=F,
      verbose=F
    );
    downregion = downregion[downregion[,"V10"]!=".",];
    downregion[,"V13"] = as.numeric(downregion[,"V13"]);
    if(sum(abs(downregion[,"V13"])<500)>0){
      downregions = downregion[abs(downregion[,"V13"])<500,c("V4","V10","V13")];
      colnames(downregions) = c("region","regionType","closestCDSDistance");
      downregions[,"regionType"] = paste(ifelse(downregions[,"closestCDSDistance"]==0,"CDS__","3UTR__"),downregions[,"regionType"],sep="");
      if(!is.null(regionType)){
        regionType = rbind(regionType,downregions);  
      } else {
        regionType = downregions;
      }
    }
  }
  sigregion.bed.merge = merge(regionType,sigregion.bed,by="region",all.y=T);
  sigregion.bed = sigregion.bed.merge[,c("chr","start","end","region","evalue","strand","gene","weightFold","regionType","closestCDSDistance")]
  #sigregion.bed[!is.na(match(siregion.bed[,"region"],regionInCDS)),""]
  #browser();
}
write.table(sigregion.bed,file="sigRegion.xls",col.names=T,row.names=F,sep="\t",quote=F);