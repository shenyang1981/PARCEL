suppressPackageStartupMessages(library("argparse"));
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

parser$add_argument("-i", "--infile", type="character", nargs=1, help="significant regions as Input", metavar="sigRegion.bed",required=T);

parser$add_argument("-g", "--genediff", type="character", nargs=1, help="differential genes", metavar="diffGenes.txt",required=T);

parser$add_argument("-s", "--siteinfo", type="character", nargs =1, help="fold change and pvalue for every site", metavar="diffSites.txt.gz",required=T);

parser$add_argument("-o", "--outdir", type="character", nargs=1, help="prefix of directory for all output", metavar="parcelResult",required=T);



parser$add_argument("-c", "--control", type="character", nargs =1, help="name of control condition in samplelist file", metavar="control",required=T);

parser$add_argument("--mincov", default=20, type="double", metavar=20,help="sites with minimum coverage of will be used in test [default %(default)s]");

parser$add_argument("--nthread", default=4, type="double", metavar=4,help="number of threads for parallel computing (multiple chromsomes) [default %(default)s]");

parser$add_argument("--batchcov", default=20, type="double", metavar=20,help="sites with minimum coverage of will be used in estimating batch effect [default %(default)s]");

parser$add_argument("--seedlogP", default=3, type="double", metavar=3,help="only sites with seed pvalue will be used for estimating Evalue [default %(default)s]");

parser$add_argument("--seedCombp", default=3, type="double", metavar=3,help="only sites with seed pvalue will be used for combined-p [default %(default)s]");

parser$add_argument("--Ecutoff", default=1, type="double", metavar=1,help="Evalue cutoff for significant region [default %(default)s]");

parser$add_argument("--topCov", default=0.01, type="double", metavar=0.01,help="mask x friction of top coverage sites [default %(default)s]");

parser$add_argument("--vsall", default="T", type="character", metavar="T",help="whether use other samples in sample list as control [default %(default)s]");

parser$add_argument("--posfile", default="", type="character", metavar="",help="bed file for positive control [default %(default)s]");

parser$add_argument("--batch", default="F", type="character", metavar="F",help="whether model batch effect [default %(default)s]");

parser$add_argument("--exactTest", default="F", type="character", metavar="F",help="whether use EdgeR exact Test [default %(default)s]");

parser$add_argument("--autobatch", default="T", type="character", metavar="T",help="whether detect batch effect automatically [default %(default)s]");

parser$add_argument("--autofilter", default="T", type="character", metavar="T",help="whether filter outlier using rosubst PCA [default %(default)s]");

parser$add_argument("--runcombp", default="T", type="character", metavar="T",help="whether run combined-p [default %(default)s]");

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
posfile = args$posfile;
cond2 = args$treatment;
cond1 = args$control;

mincov = args$mincov;
minbatchcov = args$batchcov;
seedlogP = args$seedlogP;
seedCombp = args$seedCombp;
Ecutoff = args$Ecutoff;
topCov = 1-args$topCov;
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