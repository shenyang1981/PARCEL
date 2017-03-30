suppressPackageStartupMessages(library("argparse"));
suppressPackageStartupMessages(library("bedr"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("shape"))
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


parser$add_argument("-o", "--outdir", type="character", nargs=1, help="prefix of directory for all output", metavar="parcelFilter",required=T);

parser$add_argument("-l", "--genelocusfile", type="character", nargs=1, help="exon in BED formats", metavar="NC_000964.exons.bed",required=T);

parser$add_argument("-d", "--datadir", type="character", default="", help="directory name for output files from regionEvalue.R", metavar="parcelResult");

parser$add_argument("--genediff", type="character", nargs=1, help="differential genes", metavar="diffGenes.txt",required=T);

parser$add_argument("-i", "--infile", type="character", default="", help="significant regions as Input", metavar="sigRegion.bed");
parser$add_argument("-p", "--positiveRegion", type="character", default="", help="regions of positive ressults", metavar="posRegion.bed");
parser$add_argument("-s", "--sitefile", type="character", default="", help="fold change and pvalue for every site", metavar="diffSites.txt.gz");
parser$add_argument("--Ecutoff", default=1, type="double", metavar=20,help="Evalue cutoff[default %(default)s]");
parser$add_argument("--wFCcutoff", default=0.25, type="double", metavar=20,help="weighted Fold Change cutoff[default %(default)s]");
parser$add_argument("--extWinSize", default=100, type="double", metavar=20,help="window size for calculating the wFC[default %(default)s]");
parser$add_argument("--beta", default=2, type="double", metavar=20,help="beta in cutoff function[default %(default)s]");



args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
#if ( args$verbose ) {
#  write("writing some verbose output to standard error...\n", stderr())
#}

outdir = args$outdir;
genediff = args$genediff;
datadir = args$datadir;
sigfile = args$infile;
sitefile = args$sitefile;
positivefile = args$positiveRegion;
genelocusfile = args$genelocusfile; # should be exon bed file
ecut = args$Ecutoff;
foldcut = args$wFCcutoff;
extWinSize = args$extWinSize;
beta = args$beta;

if(datadir =="" & (sitefile=="" | sigfile=="")){
  stop("datadir or significant region and site pvalues should be provided\n");
}
if(datadir!=""){
  sigfile = paste(datadir,"/","sigRegion.bed",sep="");
  sitefile = paste(datadir,"/","diffSites.txt.gz",sep="");
}

#sigfile=paste(resultdir,"sigRegion.bed",sep="/");
#sitefile=paste(resultdir,"diffSites.txt.gz",sep="/");
#"sunCount_control_sam_Batch_vsAll_10/diffSites.txt.gz";
#genediff = "~/gseq/prog/parcel/B.subtilis/samGene_control_sam_Batch_vsAll_20/geneDiff.txt";
#browser();
geneinfo= read.table(genediff,stringsAsFactors = F,quote='"',header=T);
geneinfo[,"gene"] = gsub("(.*):(.*)","\\2",geneinfo[,"id"]);
geneinfo[,"chr"] = gsub("(.*):(.*)","\\1",geneinfo[,"id"]);
#genelocusfile = "~/gseq/prog/database/Genome/UCSC/B.subtilis/NC_000964.gene.bed";
genelocus = read.table(genelocusfile,header=F,sep="\t",quote='"',stringsAsFactor=F);
genelocus = genelocus[,c(1,2,3,4,5,6)];
colnames(genelocus) = c("chr","start","end","name","score","strand");

siginfo = read.table(sigfile,stringsAsFactors = F,quote='"',header=F);
colnames(siginfo)[1:3] = c("chr","start","end");

siteinfo= read.table(sitefile,stringsAsFactors = F,quote='"',header=T);
allbin = unique(siteinfo[siteinfo[,"bin"]!="No",c("bin","Evalue")]);
siteinfo = siteinfo[,c("realchr","start","end","bin","sumcount","strand","P","logfold")];

allbin[,"chr"] = gsub("(.*):(.*):(.*)-(.*)","\\1",allbin[,"bin"]);
allbin[,"strand"] = gsub("(.*):(.*):(.*)-(.*)","\\2",allbin[,"bin"]);
allbin[allbin[,"strand"]=="Pos","strand"] = "+";
allbin[allbin[,"strand"]=="Neg","strand"] = "-";
allbin[,"start"] = as.numeric(gsub("(.*):(.*):(.*)-(.*)","\\3",allbin[,"bin"]));
allbin[,"end"] = as.numeric(gsub("(.*):(.*):(.*)-(.*)","\\4",allbin[,"bin"]));
allbin = allbin[,c("chr","start","end","bin","Evalue","strand")];
allbin[allbin[,"start"]==allbin[,"end"],"end"] = allbin[allbin[,"start"]==allbin[,"end"],"end"]+1;
allbin.sort = bedr.sort.region(allbin,check.chr=F,verbose=F,check.merge=F);
allbin.merge = bedr(
  engine = "bedtools", 
  input = list(i = allbin.sort), 
  method = "merge", 
  params = "-s -d 200 -c 5 -o min",
  check.zero.based=F,
  check.chr=F,
  verbose=F
);
colnames(allbin.merge) = c("chr","start","end","strand","Evalue");
allbin.merge[,"bin"] = bedToLocus(allbin.merge[,c("chr","start","end")]);
allbin.merge = allbin.merge[,c("chr","start","end","bin","Evalue","strand")];

positiveRegion=NULL;
if(positivefile!="" & file.exists(positivefile)){
  positiveRegion = read.table(positivefile,stringsAsFactors = F,quote='"',header=F);
  positiveRegion[,c("region")]= bedToLocus(positiveRegion[,c(1,2,3)]);
  
}

#browser();

finallistReal = getWeightFold(siteinfo = siteinfo,siginfo=siginfo,genelocus=genelocus,geneinfo=geneinfo,positiveRegion=positiveRegion,ext.win.size=extWinSize);

finallistAll = getWeightFold(siteinfo = siteinfo,siginfo=allbin.merge,genelocus=genelocus,geneinfo=geneinfo,positiveRegion=positiveRegion,ext.win.size=extWinSize)


#outdir = paste(outdir,ecut,foldcut,sep="_");
dir.create(outdir, showWarnings = F, recursive = T, mode = "0755");
setwd(outdir);

pdf("evalue_vs_wFC_all.pdf");
plot(finallistAll[,"weightFold"],finallistAll[,"logE"],pch=20,col="grey");
#mysheatmapScatter(x1 = finallistAll[,"weightFold"],x2=finallistAll[,"logE"],pch=20,ylab="-log10(Evalue)",xlab="wFC")
candpoints = finallistAll[finallistAll[,"logE"] > 1/(beta*finallistAll[,"weightFold"])^2 & finallistAll[,"maxFold"]>=1,];
if(dim(candpoints)[1]!=0){
 points(candpoints[,"weightFold"],candpoints[,"logE"],pch=20,col=adjustcolor(trop[1],alpha.f=0.9),cex=1.3);
}
if(!is.null(positiveRegion)){
 points(finallistAll[finallistAll[,"posRegion"]!=".","weightFold"],finallistAll[finallistAll[,"posRegion"]!=".","logE"],pch=20,col=adjustcolor("red",alpha.f=0.9),cex=1.1);
 
}
abline(h = -1*log10(ecut),col="grey",lty=2);
abline(v = foldcut,col="grey",lty=2);
x = seq(0,3,0.01);
y = 1/(beta*x)^2;
lines(x,y,lty=2,col="orange");
dev.off();

#browser();
#write.table(finallistReal[finallistReal[,"logE"]>=-log10(ecut) & finallistReal[,"weightFold"]>=foldcut, ],file="sigRegionFinal.xls",col.names=T,row.names=F,sep="\t",quote=F);
sigfinal = finallistReal[!is.na(finallistReal[,"weightFold"]),];
sigfinal = sigfinal[sigfinal[,"logE"] > 1/(beta*sigfinal[,"weightFold"])^2 & sigfinal[,"maxFold"]>=1, ];

write.table(sigfinal,file="sigRegionFinal.xls",col.names=T,row.names=F,sep="\t",quote=F);
write.table(finallistReal,file="sigRegionWithFC.xls",col.names=T,row.names=F,sep="\t",quote=F);
write.table(finallistAll,file="allRegionWithFC.xls",col.names=T,row.names=F,sep="\t",quote=F);
# pdf("evalue_vs_wFC_all_noBatch.pdf");
# #plot(finallistAll[,"weightFold"],finallistAll[,"logE"],pch=20,col="grey");
# mysheatmapScatter(x1 = finallistAll[,"weightFold"],x2=finallistAll[,"logE"],pch=20,ylab="-log10(Evalue)",xlab="wFC")
# points(finallistAll[finallistAll[,"posRegion"]!=".","weightFold"],finallistAll[finallistAll[,"posRegion"]!=".","logE"],pch=20,col=adjustcolor("red",alpha.f=0.9),cex=1.5);
# abline(h = -1*log10(ecut),col="grey",lty=2);
# abline(v = foldcut,col="grey",lty=2);
# dev.off();
# 
# pdf("wFC_vs_Evalue.pdf");
# plot(finallistReal[,"weightFold"],finallistReal[,"logE"]);
# if(positivefile!=""){
#   points(finallistReal[finallistReal[,"posRegion"]!=".","weightFold"],finallistReal[finallistReal[,"posRegion"]!=".","logE"],pch=20,col="red");  
# }
# dev.off();


