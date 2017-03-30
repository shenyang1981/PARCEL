suppressPackageStartupMessages(library("argparse"));
suppressPackageStartupMessages(library("bedr"))
suppressPackageStartupMessages(library("seqinr"))
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
parser$add_argument("--ismerge", default="F", type="character", metavar="F",help="whether merge libs with same condition and replicates [default %(default)s]");
parser$add_argument("--replot", default="F", type="character", metavar="F",help="only replot coverage using covinfo.xls [default %(default)s]");
parser$add_argument("--extWinSize", default=60, type="double", metavar=20,help="window size for calculating the wFC[default %(default)s]");
parser$add_argument("--genome", default="F", type="character", metavar="hg19.fas",help="genome Fasta [default %(default)s]");

parser$add_argument("--genomeSize", default="notprovided", type="character", metavar="hg19.size",help="length of chromosomes [default %(default)s]");

parser$add_argument("--getFasta", default="F", type="character", metavar="F",help="whether extract sequence from Genome [default %(default)s]");
parser$add_argument("--CDSBed", default="", type="character", metavar="CDS.bed",help="CDS annotation in BED format[default %(default)s]");
parser$add_argument("--intronBed", default="", type="character", metavar="intron.bed",help="intron annotation in BED format[default %(default)s]");
parser$add_argument("--exonBed", default="", type="character", metavar="intron.bed",help="intron annotation in BED format[default %(default)s]");

# parser$add_argument("-s", "--sitefile", type="character", default="", help="fold change and pvalue for every site", metavar="diffSites.txt.gz");
# parser$add_argument("--Ecutoff", default=1, type="double", metavar=20,help="Evalue cutoff[default %(default)s]");
# parser$add_argument("--wFCcutoff", default=0.25, type="double", metavar=20,help="weighted Fold Change cutoff[default %(default)s]");
# parser$add_argument("--extWinSize", default=100, type="double", metavar=20,help="window size for calculating the wFC[default %(default)s]");
# parser$add_argument("--beta", default=2, type="double", metavar=20,help="beta in cutoff function[default %(default)s]");



# get parameters from command line -----------------------------------


args <- parser$parse_args()

infile = args$infile;
outdir = args$outdir;
covfile = args$covfile;
samplefile = args$samplefile;
usedBatch = args$usedBatch;
cond2 = args$treatment;
cond1 = args$control;
flagVsall=as.logical(args$vsall);
ismerge = as.logical(args$ismerge);
ext.win.size = args$extWinSize;
replot = as.logical(args$replot);
getFasta = as.logical(args$getFasta);
genome = args$genome;
genomesize = args$genomeSize;
CDSBed=args$CDSBed;
intronBed=args$intronBed;
exonBed=args$exonBed;
# infile="testevalueFilter_1_0.25/sigRegionFinal.xls"
# covfile="testevalue_batch1_control_sam_Batch_vsAll_20/all.cov.txt"
# samplefile = "sampleList_Bsub.txt";
# usedBatch = "batch1";
# cond2 = "sam";
# cond1 = "control";
# outdir = "testout";
# flagVsall=T;
# ext.win.size=100;
#browser();

# load sample information -------------------------------------------------

chrsize = read.table(genomesize,header=F,sep="\t",quote='"');
colnames(chrsize) = c("chr","size");

sampleInfo = read.table(samplefile,header=T,sep="\t",stringsAsFactors = F);

sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == usedBatch,];

if(ismerge){
  sampleInfo[,"LibID"] = paste(sampleInfo[,"Condition"],sampleInfo[,"Replicates"],sep="__");
}
tmpSampleInfo = unique(sampleInfo[,c("Condition","LibID")]);

if(flagVsall==T){
  tmpSampleInfo[tmpSampleInfo[,"Condition"]!=cond2,"Condition"] = cond1;  
}

treatmentCols = tmpSampleInfo[tmpSampleInfo[,"Condition"]==cond2,"LibID"];
controlCols = tmpSampleInfo[tmpSampleInfo[,"Condition"]==cond1,"LibID"];


# get CDS regions from bed file -------------------------------------------
#browser();
sigregion = read.table(infile,header=T,sep="\t",stringsAsFactors = F)

# only keep chrs from sig regions
usedchrs = unique(sigregion[,"chr"]);

cds.sort = NULL;
if(file.exists(CDSBed)){
  cds = read.table(CDSBed,header=F,sep="\t",quote='"',stringsAsFactor=F);
  cds = cds[!is.na(match(cds[,1],usedchrs)),];
  
  cds.sort = bedr.sort.region(cds,check.chr=F,check.zero.based = F,verbose = F);
  
  first.cds = getEntryFromOrderBed(x = cds.sort,myorder="first");
  
  first.cds.sort = bedr.sort.region(first.cds,check.chr=F,check.zero.based = F,verbose = F);
  
  last.cds = getEntryFromOrderBed(x = cds.sort,myorder="last");
  
  last.cds.sort = bedr.sort.region(last.cds,check.chr=F,check.zero.based = F,verbose = F);
  
  cds.length = aggregate(cds.sort[,3]-cds.sort[,2]+1,list(cds.sort[,4]),sum);
  
  colnames(cds.length) = c("cdsgene","CDSLength");
}

dir.create(outdir, showWarnings = F, recursive = T, mode = "0755");

if(dim(sigregion)[1]==0){
  setwd(outdir);
  # has to output something in case the snakemake is down.
  writeLines("", con = "sequences.fas");
  writeLines("", con = "sequencesRegionID.clean.fas");
  writeLines("", con = "sigRegion.xls");
  writeLines("", con = "covinfo.xls");
  cat("No Candidates\n");
  quit("no");
}

if(replot==F){
  sigregion.bed = sigregion[,c("chr","start","end","region","evalue","strand","gene","weightFold")]
  sigregion.bed[,"region"] = bedToLocus(sigregion.bed[,c("chr","start","end")]);
  sigregion.bed = merge(sigregion.bed,chrsize,by="chr");
  
  sigMidpoint = round((sigregion.bed[,2] + sigregion.bed[,3])/2,0);
  sigWidth = sigregion.bed[,3] - sigregion.bed[,2] + 1;
  if(sum(sigWidth<ext.win.size*2)>0){
    sigregion.bed[sigWidth<ext.win.size*2,2] = ifelse(sigMidpoint[sigWidth<ext.win.size*2] -ext.win.size>0,sigMidpoint[sigWidth<ext.win.size*2] -ext.win.size,1);
    sigregion.bed[sigWidth<ext.win.size*2,3] = sigMidpoint[sigWidth<ext.win.size*2] +ext.win.size;
  }
  
  sigregion.bed[sigregion.bed[,"start"]<1,"start"] = 1;
  sigregion.bed[sigregion.bed[,"end"]>sigregion.bed[,"size"],"end"] = sigregion.bed[sigregion.bed[,"end"]>sigregion.bed[,"size"],"size"];
  
  sigregion.sort = bedr.sort.region(sigregion.bed,check.chr=F,check.zero.based = F,verbose = F)
  
  sites = read.table(covfile,header=T,sep="\t",quote='"',comment.char='$',stringsAsFactors=F)
  # tmp debug should be deleted later
  tmp1 = colnames(sites);
  tmp1 = gsub("\\.","-",tmp1);
  colnames(sites) = tmp1;
  
  #
  sites = sites[!is.na(match(sites[,1],usedchrs)),];
  
  sites.sort = bedr.sort.region(sites,check.chr=F,check.zero.based = F,verbose = F,check.merge = F)
  
  #sigsites = bedr.join.region(sigregion.sort,sites.sort,check.chr=F,check.merge = F,check.zero.based = F,verbose=F)
  #browser();
  sigsites = bedr(
    engine = "bedtools", 
    input = list(a = sigregion.sort ,b=sites.sort), 
    method = "intersect", 
    params = "-s -loj -sorted",
    check.chr=F,
    verbose=F
  );
  sigsites = sigsites[sigsites[,"#chr.b"]!=".",];
  if(dim(sigsites)[1]==0){
    cat("Error, no sites coverage\n");
    quit("no");
  }
  if(!is.null(cds.sort)){
    tmpsig = sigsites[,c("#chr.b","start.b","end.b","name","score","strand")];
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
    sigsitesIndex = bedToLocus(sigsites[,c("#chr.b","start.b","end.b")]);
    cdssiteIndex = bedToLocus(tmpCDSsites[tmpCDSsites[,"V4"]!=".",c("#chr.b","start.b","end.b")]);
    sigsites[,"isCDS"] = "No";
    sigsites[!is.na(match(sigsitesIndex,cdssiteIndex)),"isCDS"] = "Yes";
  }
  #browser();
  sigsites[,"evalue"] = as.numeric(sigsites[,"evalue"]);
  sigsites[,"weightFold"] = as.numeric(sigsites[,"weightFold"]);
  sigsites[,"minchange"] = as.numeric(sigsites[,"minchange"]);
  
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
  sigsites = sigsites[!is.na(match(sigsites[,"chr"],usedchrs)),];
  
  sigregion.bed = sigregion[,c("chr","start","end","region","evalue","strand","gene","weightFold")]
  sigregion.bed[,"region"] = bedToLocus(sigregion.bed[,c("chr","start","end")]);
  sigregion.bed = merge(sigregion.bed,chrsize,by="chr");
  
  sigMidpoint = round((sigregion.bed[,2] + sigregion.bed[,3])/2,0);
  sigWidth = sigregion.bed[,3] - sigregion.bed[,2] + 1;
  
  if(sum(sigWidth<ext.win.size*2)>0){
    sigregion.bed[sigWidth<ext.win.size*2,2] = sigMidpoint[sigWidth<ext.win.size*2] -ext.win.size;
    sigregion.bed[sigWidth<ext.win.size*2,3] = sigMidpoint[sigWidth<ext.win.size*2] +ext.win.size;
  }
  
  sigregion.bed[sigregion.bed[,"start"]<1,"start"] = 1;
  sigregion.bed[sigregion.bed[,"end"]>sigregion.bed[,"size"],"end"] = sigregion.bed[sigregion.bed[,"end"]>sigregion.bed[,"size"],"size"];
  
  sigregion.bed = sigregion.bed[,c("chr","start","end","region","evalue","strand","gene","weightFold","region")]
}


# get intron region from bed file -----------------------------------------

intron.sort = NULL;
if(file.exists(intronBed)){
  intron = read.table(intronBed,header=F,sep="\t",quote='"',stringsAsFactor=F);
  intron = intron[!is.na(match(intron[,1],usedchrs)),];
  
  intron.sort = bedr.sort.region(intron,check.chr=F,check.zero.based = F,verbose = F);
}

# get exon region from bed file -----------------------------------------

exon.sort = NULL;
if(file.exists(exonBed)){
  exon = read.table(exonBed,header=F,sep="\t",quote='"',stringsAsFactor=F);
  exon = exon[!is.na(match(exon[,1],usedchrs)),];
  
  exon.sort = bedr.sort.region(exon,check.chr=F,check.zero.based = F,verbose = F);
}

# set output dir ----------------------------------------------------------

dir.create(outdir, showWarnings = F, recursive = T, mode = "0755");


# load candidate regions --------------------------------------------------

# set outdir
#outdir = paste(outdir,ecut,foldcut,sep="_");
#browser();
#ddply(sigsites,.(region),plotCoverage,treatmentCols=treatmentCols,controlCols=controlCols,genome=genome,getFasta=getFasta,cds.sort=cds.sort);
#browser();

fillzero = function(x=NULL,fillcols=NULL){
  # fill zero to those sites with coverage 0
  tmpx = seq(x[1,"start"]-1,x[1,"end"]+1,1);
  tmpnofillcols = x[,is.na(match(colnames(x),fillcols))];
  rangestart = as.numeric(gsub("(.*):(.*)-(.*)","\\2",x[1,"region"]));
  rangeend = as.numeric(gsub("(.*):(.*)-(.*)","\\3",x[1,"region"]));
  if(x[1,"strand"]=="+"){
    relativePos = rangestart - tmpx + 1;
  } else {
    relativePos = rangeend - tmpx + 1;
  }
  #tmpdata = data.frame("start.b"=tmpx);
  #all = merge(tmpnofillcols,tmpdata,by="start.b",all.y=T);
  tmpdata = data.frame("chr" = x[1,"chr"],"start" = x[1,"start"]-1,"end"=x[1,"end"]+1,"region"=x[1,"region"],"evalue"=x[1,"evalue"],"strand"=x[1,"strand"],"gene"=x[1,"gene"],"start.b"=tmpx,"relativePos"=relativePos,"weightFold"=0,"score"=0,"minchange"=0,stringsAsFactors = F);
  rownames(tmpdata) = as.character(tmpx);
  tmpdata[,fillcols] = 0;
  tmpdata[as.character(x[,"start.b"]),c("weightFold","score","minchange",fillcols)] = x[,c("weightFold","score","minchange",fillcols)];
  if(sum(is.na(tmpdata[,"region"]))>0){
    browser();
  }
  return(tmpdata);
}
#browser();
covinfo = sigsites[,c("chr","start","end","region","evalue","strand","gene","weightFold","start.b","score","minchange",treatmentCols,controlCols)];
covinfo = ddply(covinfo,.(region),fillzero,fillcols = c(treatmentCols,controlCols));
covinfo = covinfo[order(covinfo[,"region"],covinfo[,"relativePos"]),];

# get the minchange for region
regionminchange = function(x){
  result = x[order(x[,"score"],decreasing = T),][1,c("minchange","start.b")];
  return(result);
}
minchange = ddply(covinfo,.(region),regionminchange);
colnames(minchange) = c("region","minchange","sigsite");

# get fasta sequence
#browser();
regionLowComplex = NULL;
if(getFasta==T & file.exists(genome) & file.exists(genomesize)){
  # get sequence from region block,filter low complexity sequence
  tmpcovRegion = unique(covinfo[,c("chr","start","end","region","evalue","strand")]);
  tmpcovRegion[,c("chr","start","end")] = index2bed(tmpcovRegion[,"region"]);
  tmpMidpoint = round((tmpcovRegion[,2] + tmpcovRegion[,3])/2,0);
  tmp.ext.win.size = 24;
  tmpcovRegion[,3] = ifelse(tmpcovRegion[,"strand"]=="+", tmpMidpoint + tmp.ext.win.size*2,tmpMidpoint);
  tmpcovRegion[,2] = ifelse(tmpcovRegion[,"strand"]=="-", tmpMidpoint - tmp.ext.win.size*2,tmpMidpoint);
  
  tmpcovRegion = merge(tmpcovRegion,chrsize,by="chr");
  
  tmpcovRegion[tmpcovRegion[,"start"]<1,"start"] = 1;
  tmpcovRegion[tmpcovRegion[,"end"]>tmpcovRegion[,"size"],"end"] = tmpcovRegion[tmpcovRegion[,"end"]>tmpcovRegion[,"size"],"size"];
  
  tmpseq = get.fasta(tmpcovRegion, fasta=genome, output.fasta = F,check.chr=F,check.sort = F,check.merge = F,check.zero.based = F,strand=T,verbose=F,use.name.field = T);
  tmpseq[,"Lcomplex"] = unlist(lapply(tmpseq[,2],function(x){
    nucCount = table(unlist(strsplit(toupper(x),split="")));
    if(max(nucCount/sum(nucCount))>0.5){
      return(T);
    } else {
      return(F);
    }
  }));
  regionLowComplex = tmpseq[,c("index","Lcomplex")];
  colnames(regionLowComplex) = c("region","Lcomplex");
  
  # write sequence without filtering
  tmpcovRegion = unique(covinfo[,c("chr","start","end","region","evalue","strand")]);
  
  tmpcovRegion = merge(tmpcovRegion,chrsize,by="chr");
  tmpcovRegion[tmpcovRegion[,"start"]<1,"start"] = 1;
  tmpcovRegion[tmpcovRegion[,"end"]>tmpcovRegion[,"size"],"end"] = tmpcovRegion[tmpcovRegion[,"end"]>tmpcovRegion[,"size"],"size"];
  
  tmpseq = get.fasta(tmpcovRegion, fasta=genome, output.fasta = TRUE,check.chr=F,check.sort = F,check.merge = F,check.zero.based = F,strand=T,verbose=F);
  writeLines(tmpseq[[1]], con = "sequences.fas")
  
  # write final sequence
  tmpcovRegion = unique(covinfo[,c("chr","start","end","region","evalue","strand")]);
  tmpcovRegion = tmpcovRegion[is.na(match(tmpcovRegion[,"region"],regionLowComplex[regionLowComplex[,"Lcomplex"]==T,"region"])),];
  if(dim(tmpcovRegion)[1]!=0){
    tmpcovRegion = merge(tmpcovRegion,chrsize,by="chr");
    tmpcovRegion[tmpcovRegion[,"start"]<1,"start"] = 1;
    tmpcovRegion[tmpcovRegion[,"end"]>tmpcovRegion[,"size"],"end"] = tmpcovRegion[tmpcovRegion[,"end"]>tmpcovRegion[,"size"],"size"];
    
    tmpseq = get.fasta(tmpcovRegion, fasta=genome, output.fasta = TRUE,check.chr=F,check.sort = F,check.merge = F,check.zero.based = F,strand=T,verbose=F);
    writeLines(tmpseq[[1]], con = "sequences.clean.fas")
    
    tmpseqid = get.fasta(bedr.sort.region(tmpcovRegion,check.chr=F,check.zero.based = F,verbose = F,check.merge=F), fasta=genome, output.fasta = TRUE,check.chr=F,check.sort = F,check.merge = F,check.zero.based = F,strand=T,verbose=F,use.name.field = T);
    writeLines(tmpseqid[[1]], con = "sequencesRegionID.clean.fas")  
    
    tmpcovRegion[,"region"] = paste(tmpcovRegion[,"region"],paste(bedToLocus(tmpcovRegion[,1:3]),"(",tmpcovRegion[,"strand"],")",sep=""),sep="__");
    tmpseqid = get.fasta(bedr.sort.region(tmpcovRegion,check.chr=F,check.zero.based = F,verbose = F,check.merge=F), fasta=genome, output.fasta = TRUE,check.chr=F,check.sort = F,check.merge = F,check.zero.based = F,strand=T,verbose=F,use.name.field = T);
    writeLines(tmpseqid[[1]], con = "sequencesRegionIDLocus.clean.fas")  
  }
}
write.table(covinfo,file="covinfo.xls",col.names=T,row.names=F,sep="\t",quote=F);


# annotate biotype for different region -----------------------------------

if(!is.null(cds.sort)){
  sigMidpoint = round((sigregion.sort[,2] + sigregion.sort[,3])/2,0);
  sigregion.sort[,2] =  sigMidpoint;
  sigregion.sort[,3] =  sigMidpoint+1;
  #browser();
  #regionType = NULL;
  
  #first see whether it located in a upstream region of any first exon of CDS  within 500bp
  upregion = bedr(
      engine = "bedtools", 
      input = list(a = sigregion.sort[,1:6] ,b=first.cds.sort[,1:6]), 
      method = "closest", 
      params = "-s -D b -d -io -id",
      check.chr=F,
      verbose=F
    );
  upregion[,"V13"] = as.numeric(upregion[,"V13"]);
  upregion = upregion[upregion[,"V10"]!="." & abs(upregion[,"V13"])<500,];
  regionanno=NULL;
  if(dim(upregion)[1]>0){
    regionanno = upregion[,c("V4","V10","V13")];
    colnames(regionanno) = c("region","regionType","closestCDSDistance");
    regionanno[,"regionType"] = paste("5UTR__",regionanno[,"regionType"],sep="");
  }
  
  #sigregion.sort.new = sigregion.sort[is.na(match(sigregion.sort[,"region"],regionanno[,"region"])),];
  #browser();
  innerregion = bedr(
    engine = "bedtools", 
    input = list(a = sigregion.sort[,1:6],b=cds.sort[,1:6]), 
    method = "intersect", 
    params = "-s -loj -sorted",
    check.chr=F,
    verbose=F
  );
  if(!is.null(regionanno)){
    innerregion = innerregion[innerregion[,"V4"]!="." & is.na(match(innerregion[,"region"],regionanno[,"region"])),];  
  } else {
    innerregion = innerregion[innerregion[,"V4"]!=".",];
  }
  
  if(dim(innerregion)[1]!=0){
    innerregion[,c("V2","V3","V5","start","end")] = apply(innerregion[,c("V2","V3","V5","start","end")],2,as.numeric);
    tmpdistance = ifelse(innerregion[,"V6"]=="+",innerregion[,"start"]-innerregion[,"V2"]+innerregion[,"V5"],innerregion[,"V3"]-innerregion[,"end"]+innerregion[,"V5"]);
    tmp = data.frame("region" = innerregion[,"region"],"regionType"=innerregion[,"V4"],"closestCDSDistance"=tmpdistance,stringsAsFactors = F);
    tmp[,"regionType"] = paste("CDS__",tmp[,"regionType"],sep="");
    regionanno = rbind(regionanno,tmp);
  }
  
  downregion = bedr(
    engine = "bedtools", 
    input = list(a = sigregion.sort[,1:6] ,b=last.cds.sort[,1:6]), 
    method = "closest", 
    params = "-s -D b -d -io -iu",
    check.chr=F,
    verbose=F
  );
  downregion[,"V13"] = as.numeric(downregion[,"V13"]);
  if(!is.null(regionanno)){
    downregion = downregion[downregion[,"V10"]!="." & abs(downregion[,"V13"])<500 & is.na(match(downregion[,"V4"],regionanno[,"region"])),];  
  } else {
    downregion = downregion[downregion[,"V10"]!=".",];
  }
  
  if(dim(downregion)[1]>0){
    tmp = downregion[,c("V4","V10","V13")];
    colnames(tmp) = c("region","regionType","closestCDSDistance");
    tmp[,"regionType"] = paste("3UTR__",tmp[,"regionType"],sep="");
    regionanno = rbind(regionanno,tmp);
  }

  # annotate exonic region
  if(!is.null(exon.sort)){
    exonregion = bedr(
      engine = "bedtools", 
      input = list(a = sigregion.sort[,1:6],b=exon.sort[,1:6]), 
      method = "intersect", 
      params = "-s -loj -sorted",
      check.chr=F,
      verbose=F
    );
    if(!is.null(regionanno)){
      exonregion = exonregion[exonregion[,"V4"]!="." & is.na(match(exonregion[,"region"],regionanno[,"region"])),];  
    } else {
      exonregion = exonregion[exonregion[,"V4"]!=".",];
    }
    
    if(dim(exonregion)[1]!=0){
      exonregion[,c("V2","V3","V5","start","end")] = apply(exonregion[,c("V2","V3","V5","start","end")],2,as.numeric);
      tmp = data.frame("region" = exonregion[,"region"],"regionType"=exonregion[,"V4"],"closestCDSDistance"=0,stringsAsFactors = F);
      tmp[,"regionType"] = paste("Exon__",tmp[,"regionType"],sep="");
      regionanno = rbind(regionanno,tmp);
    }
  }
  
  # annotate intronic region
  if(!is.null(intron.sort)){
    intronregion = bedr(
      engine = "bedtools", 
      input = list(a = sigregion.sort[,1:6],b=intron.sort[,1:6]), 
      method = "intersect", 
      params = "-s -loj -sorted",
      check.chr=F,
      verbose=F
    );
    if(!is.null(regionanno)){
      intronregion = intronregion[intronregion[,"V4"]!="." & is.na(match(intronregion[,"region"],regionanno[,"region"])),];  
    }
    if(dim(intronregion)[1]!=0){
      intronregion[,c("V2","V3","V5","start","end")] = apply(intronregion[,c("V2","V3","V5","start","end")],2,as.numeric);
      tmp = data.frame("region" = intronregion[,"region"],"regionType"=intronregion[,"V4"],"closestCDSDistance"=0,stringsAsFactors = F);
      tmp[,"regionType"] = paste("Intron__",tmp[,"regionType"],sep="");
      regionanno = rbind(regionanno,tmp);
    }
  }
    
  if(is.null(regionanno)){
    regionanno = data.frame("region"="","regionType"="","closestCDSDistance"=0);
  }
  sigregion.bed.merge = merge(regionanno,sigregion.bed,by="region",all.y=T);
  sigregion.bed.merge[,"cdsgene"] = gsub("(.*)__(.*)","\\2",sigregion.bed.merge[,"regionType"]);
  sigregion.bed.merge = merge(sigregion.bed.merge,cds.length,by="cdsgene",all.x=T);
  sigregion.bed = sigregion.bed.merge[,c("chr","start","end","region","evalue","strand","gene","weightFold","regionType","closestCDSDistance","CDSLength")]
  
}

# loading sequence
seqs = "sequencesRegionIDLocus.clean.fas";
myseqs = read.fasta(seqs);
myseqtab = data.frame("seqs" = unlist(lapply(myseqs,function(x){toupper(c2s(getSequence(x)))})),"region" = unlist(lapply(myseqs,function(x){getName(x)})));
myseqtab[,"seqLocus"] = gsub(".*__(.*)","\\1",myseqtab[,"region"]);
myseqtab[,"region"] = gsub("(.*)__(.*)","\\1",myseqtab[,"region"]);
sigregion.bed = merge(sigregion.bed,myseqtab,by="region",all.x=T);
sigregion.bed = merge(sigregion.bed,minchange,by="region");
sigregion.bed = merge(sigregion.bed,regionLowComplex,by="region");
write.table(sigregion.bed,file="sigRegion.xls",col.names=T,row.names=F,sep="\t",quote=F);
write.table(sigregion.bed[,1:6],file="sigRegion.bed",col.names=F,row.names=F,sep="\t",quote=F);
