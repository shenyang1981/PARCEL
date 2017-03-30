# compare results between pipeline
options(stringsAsFactors=F);
library("VennDiagram");
library("bedr");
#setwd("~/gseq/prog/parcel/B.subtilis/")

args = commandArgs(T);
infile1 = args[1];
infile2 = args[2];
output = args[3];
#a = read.table("sumResult/Bsub_SAM_v2all_output2_wFilters.bed",header=T,sep="\t",quote='"');
#b = read.table("samFilter_NoMask_batch1_1_0.25/sigRegionFinal.xls",header=T,sep="\t",quote='"');

a = read.table(infile1,header=T,sep="\t",quote='"');
b = read.table(infile2,header=T,sep="\t",quote='"')
a[,"region"] = paste(a[,"region"],"alg1",sep="__");
b[,"region"] = paste(b[,"region"],"alg2",sep="__");

all = rbind(a[,c("chr","start","end","region","evalue","strand")],b[,c("chr","start","end","region","evalue","strand")])

all.sort = bedr.sort.region(all,check.chr=F,check.zero.based = F,check.valid = F,check.merge = F);

all.merge = bedr(
  engine = "bedtools", 
  input = list(i = all.sort), 
  method = "merge", 
  params = "-s -d 200 -c 4,5 -o distinct,min",
  check.zero.based=F,
  check.chr=F,
  verbose=F
);

colnames(all.merge) = c("chr","start","end","strand","region","evalue");
total = length(unique(all.merge[,"region"]));
regalg1 = length(unique(all.merge[grepl("alg1",all.merge[,"region"]) & !grepl("alg2",all.merge[,"region"]),"region"]));
regalg2 = length(unique(all.merge[grepl("alg2",all.merge[,"region"]) & !grepl("alg1",all.merge[,"region"]),"region"]));
regall = length(unique(all.merge[grepl("alg2",all.merge[,"region"]) & grepl("alg1",all.merge[,"region"]),"region"]));

#plot(venneuler(A=regalg1,B=regalg2,"A&B"=regall),main=paste("Total : ", total,"\n","Shared: ",regall,sep=""));
pdf(output);
draw.pairwise.venn(regalg1+regall, regalg2+regall, regall, category = c("Alg1", "Alg2"), 
                   lty = rep("blank",2), fill = c("light blue", "pink"), cex=4,
                   alpha = rep(0.5, 2), cat.pos = c(0,0), cat.dist = rep(0.025, 2));
dev.off();
