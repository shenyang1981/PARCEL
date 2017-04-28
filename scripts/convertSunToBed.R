suppressPackageStartupMessages(library("argparse"));
suppressPackageStartupMessages(library("bedr"));

parser <- ArgumentParser(description='Plot the coverage for significant region');

parser$add_argument("-i", "--infile", type="character", nargs=1, help="candidates from SunMiao's Pipeline as Input", metavar="combined_lys_output2_wfilters.txt",required=T);

parser$add_argument("-s", "--candseqs", type="character", nargs=1, help="output candidates sequence to here", metavar="candseqs.fa",required=T);

parser$add_argument("-o", "--output", type="character", nargs=1, help="converted region with genomic coordinates", metavar="combined_met_output2_wfilters.bed",required=T);

parser$add_argument("-r", "--ref", type="character", nargs =1, help="reference transcriptome files", metavar="sgdGenes_addUTR_spiked1.fa",required=T);
parser$add_argument("-g", "--genome", type="character", nargs =1, help="reference genome files", metavar="genome.fa",required=T);

#parser$add_argument("--mincov", default=20, type="double", metavar=20,help="sites with minimum coverage of will be used in test [default %(default)s]");

#parser$add_argument("--nthread", default=4, type="double", metavar=4,help="number of threads for parallel computing (multiple chromsomes) [default %(default)s]");

args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
#if ( args$verbose ) {
#  write("writing some verbose output to standard error...\n", stderr())
#}


suncandfile = args$infile;
refTt = args$ref;
refGenome = args$genome;
candseqs = args$candseqs;
output = args$output;

sunResult=read.table(suncandfile,header=T,sep="\t",quote='"');
if(dim(sunResult)[1]==0){
  #result = data.frame("chr"=NULL,"start"=NULL,"end"=NULL,"region"=NULL,"evalue"=NULL,"strand"=NULL,"gene"=NULL,"weightFold"=NULL);
  result = data.frame("chr","start","end","region","evalue","strand","gene","weightFold");
  write.table(result,file=output,col.names=F,row.names=F,sep="\t",quote=F);
  quit("no");
}
sunbed = sunResult[,c("geneID","pstart","pend","E_value","maxFC")];
sunbed[,"region"] = paste(sunbed[,"geneID"],":",sunbed[,"pstart"],"-",sunbed[,"pend"],sep="");
sunbed[,"strand"] = "+";
sunbed[,"chr"] = sunbed[,"geneID"];
sunbed[,"gene"] = sunbed[,"geneID"];
sunbed[,"weightFold"] = sunbed[,"maxFC"];
sunbed[,"evalue"] = sunbed[,"E_value"];
sunbed[,"start"] = sunbed[,"pstart"];
sunbed[,"end"] = sunbed[,"pend"];

sunbed = sunbed[,c("chr","start","end","region","evalue","strand","gene","weightFold")];
sunbed.sort = bedr.sort.region(sunbed,check.chr=F,check.zero.based = F,check.valid = F,check.merge = F,engine = "bedtools");
sunbed.sort[,"oldstart"] = sunbed.sort[,"start"];
sunbed.sort[,"oldend"] = sunbed.sort[,"end"];
sunbed.sort[,"start"] = ifelse(sunbed.sort[,"start"]>50,sunbed.sort[,"start"]-50,sunbed.sort[,"start"]);
sunbed.sort[,"end"] = ifelse(sunbed.sort[,"start"]<50,sunbed.sort[,"end"]+50,sunbed.sort[,"end"]);
sunbed.sort[,"region"] = paste(sunbed.sort[,"gene"],":",sunbed.sort[,"start"],"-",sunbed.sort[,"end"],sep="");

sunseqs = get.fasta(sunbed.sort,fasta = refTt,strand = T,check.zero.based = F,check.chr=F,check.valid=F,check.merge=F,output.fasta = T)
writeLines(sunseqs[[1]], con = candseqs)

blatcmd = paste("blat",refGenome,candseqs,"stdout | perl ~/gseq/prog/script/psl_to_bed_best_score.pl 0.95",sep=" ");
blatresult = system(blatcmd,intern = T);

if(length(blatresult)==0){
  result = data.frame("chr","start","end","region","evalue","strand","gene","weightFold");
  write.table(result,file=output,col.names=F,row.names=F,sep="\t",quote=F);
  quit("no");
}

blattab = t(matrix(sapply(blatresult,function(x){unlist(strsplit(x,split="\t"))}),nrow=13));
if(dim(blattab)[1]==1){
  blattab = data.frame(t(blattab[,c(1,4,7,2,3,8,9)]));
} else {
  blattab = as.data.frame(blattab[,c(1,4,7,2,3,8,9)]);
}
blattab[,2] = gsub("(>*)\\(.*","\\1",blattab[,2]);
colnames(blattab) = c("chr","region","gstrand","gstart","gend","tstart","tend");
blattab[,c(-1,-2,-3)] = apply(blattab[,c(-1,-2,-3)],2,function(x){as.numeric(as.character(x))});
sunbed.sort.new = merge(sunbed.sort,blattab,by="region")

# convert the coordinate to locus before extension
sunbed.sort.new[,"blocksize"] = (sunbed.sort.new[,"oldend"] - sunbed.sort.new[,"oldstart"] + 1); # only if the original block is located inside the aligned block
sunbed.sort.new[,"startshift"] = (sunbed.sort.new[,"oldstart"] - sunbed.sort.new[,"start"]);
# pos strand
tmpindex = sunbed.sort.new[,"gstrand"]=="+";
if(sum(tmpindex)>0){
  sunbed.sort.new[tmpindex,"ogstart"] = sunbed.sort.new[tmpindex,"gstart"] + sunbed.sort.new[tmpindex,"startshift"] - sunbed.sort.new[tmpindex,"tstart"];
  sunbed.sort.new[tmpindex,"ogend"] = sunbed.sort.new[tmpindex,"ogstart"] + sunbed.sort.new[tmpindex,"blocksize"] - 1;
}
tmpindex = sunbed.sort.new[,"gstrand"]=="-";
if(sum(tmpindex)>0){
  sunbed.sort.new[tmpindex,"ogend"] = sunbed.sort.new[tmpindex,"gend"] - sunbed.sort.new[tmpindex,"startshift"] + sunbed.sort.new[tmpindex,"tstart"];
  sunbed.sort.new[tmpindex,"ogstart"] = sunbed.sort.new[tmpindex,"ogend"] - sunbed.sort.new[tmpindex,"blocksize"] + 1;
}

sunbed.sort.new[,"region"] = paste(sunbed.sort.new[,"chr.y"],":",sunbed.sort.new[,"ogstart"],"-",sunbed.sort.new[,"ogend"],sep="");
sunbed.sort.out = sunbed.sort.new[,c("chr.y","ogstart","ogend","region","evalue","gstrand","gene","weightFold")];
colnames(sunbed.sort.out) = c("chr","start","end","region","evalue","strand","gene","weightFold");
write.table(sunbed.sort.out,file=output,col.names=T,row.names=F,sep="\t",quote=F);
