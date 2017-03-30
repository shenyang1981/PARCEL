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

parser$add_argument("-g", "--genediff", type="character", nargs=1, help="differential genes", metavar="diffGenes.txt",required=T);

parser$add_argument("-l", "--genelocusfile", type="character", nargs=1, help="exon in BED formats", metavar="NC_000964.exons.bed",required=T);

parser$add_argument("-d", "--datadir", type="character", nargs=1, help="directory name for output files from regionEvalue.R", metavar="parcelResult",required=F);

parser$add_argument("-i", "--infile", type="character", nargs=1, help="significant regions as Input", metavar="sigRegion.bed",required=F);

parser$add_argument("-p", "--positiveRegion", type="character", nargs=1, help="regions of positive ressults", metavar="posRegion.bed",required=F);

parser$add_argument("-s", "--sitefile", type="character", nargs =1, help="fold change and pvalue for every site", metavar="diffSites.txt.gz",required=F);




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


sigfile=paste(resultdir,"sigRegion.bed",sep="/");
sitefile=paste(resultdir,"diffSites.txt.gz",sep="/");
#"sunCount_control_sam_Batch_vsAll_10/diffSites.txt.gz";
genediff = "~/gseq/prog/parcel/B.subtilis/samGene_control_sam_Batch_vsAll_20/geneDiff.txt";

geneinfo= read.table(genediff,stringsAsFactors = F,quote='"',header=T);
geneinfo[,"gene"] = gsub("(.*):(.*)","\\2",geneinfo[,"id"]);

genelocusfile = "~/gseq/prog/database/Genome/UCSC/B.subtilis/NC_000964.gene.bed";
genelocus = read.table(genelocusfile,header=F,sep="\t",quote='"',stringsAsFactor=F);
genelocus = genelocus[,c(1,2,3,4,5,6)];
colnames(genelocus) = c("chr","start","end","name","score","strand");


siginfo = read.table(sigfile,stringsAsFactors = F,quote='"',header=F);
colnames(siginfo)[1:3] = c("chr","start","end");

siteinfo= read.table(sitefile,stringsAsFactors = F,quote='"',header=T);
allbin = unique(siteinfo[siteinfo[,"bin"]!="No",c("bin","Evalue")]);
siteinfo = siteinfo[,c("chr","start","end","P","logfold","sumcount")];

allbin[,"chr"] = gsub("(.*):(.*)-(.*)","\\1",allbin[,"bin"]);
allbin[,"start"] = as.numeric(gsub("(.*):(.*)-(.*)","\\2",allbin[,"bin"]));
allbin[,"end"] = as.numeric(gsub("(.*):(.*)-(.*)","\\3",allbin[,"bin"]));
allbin = allbin[,c("chr","start","end","bin","Evalue")];
allbin[allbin[,"start"]==allbin[,"end"],"end"] = allbin[allbin[,"start"]==allbin[,"end"],"end"]+1;
allbin.sort = bedr.sort.region(allbin,check.chr=F,verbose=F);
allbin.merge = bedr(
  engine = "bedtools", 
  input = list(i = allbin.sort), 
  method = "merge", 
  params = "-d 200 -c 5 -o min",
  check.zero.based=F,
  check.chr=F,
  verbose=F
);
colnames(allbin.merge) = c("chr","start","end","Evalue");
allbin.merge[,"bin"] = bedToLocus(allbin.merge[,c("chr","start","end")]);
allbin.merge = allbin.merge[,c("chr","start","end","bin","Evalue")];

positiveRegion = read.table(positivefile,stringsAsFactors = F,quote='"',header=F);
positiveRegion[,c("region")]= bedToLocus(positiveRegion[,c(1,2,3)]);

