#setwd("~/gseq/prog/parcel/B.subtilis/")
suppressPackageStartupMessages(library("bedr"))

args = commandArgs(T);
infile = args[1];
outfile = args[2];
ismerge = args[3];
posfile = args[4];

# infile = "~/gseq/prog/parcel/B.subtilis/sumResult/Bsub_SAM_v2all_output2_wFilters.txt"
# outfile = "tmp.txt"
# posfile = "~/gseq/prog/parcel/RegPrecise/Bsub__SAM.bed"

a = read.table(infile,header=T,sep="\t",quote='"',stringsAsFactors=F)
b = a[,c("pos","strand","E_value","maxFC","geneID","winSize")];
b[,"chr"] = "chr";
b[,"start"] = b[,"pos"]-round(b[,"winSize"]/2,0);
b[,"end"] = b[,"pos"] + b[,"winSize"];
b[,"region"] = paste(b[,"chr"],":",b[,"start"],"-",b[,"end"],sep="");
b[,"weightFold"] = b[,"maxFC"];
b[,"gene"] = b[,"geneID"];
b[,"evalue"] = b[,"E_value"];
b = b[,c("chr","start","end","region","evalue","strand","gene","weightFold")];
result = b;
if(file.exists(posfile)){
  posbed = read.table(posfile,header=F,sep="\t",quote='"',stringsAsFactors=F);
  colnames(posbed) = c("chr","start","end","region","evalue","strand");
  posbed[,"gene"] = gsub("(.*)\\(.*","\\1",posbed[,"region"]);
  posbed[,"gene"] = paste(posbed[,"gene"],"__Positive",sep="");
  posbed[,"weightFold"] = 0;
  posbed[,"evalue"]=100;
  allbed = rbind(b,posbed);
  allbed.sort = bedr.sort.region(allbed,check.chr = F,check.valid = F,check.merge = F,verbose=F);
  if(!is.null(ismerge) & ismerge == "merge"){
    allbed.sort.merge = bedr(
      engine = "bedtools", 
      input = list(i = allbed.sort), 
      method = "merge", 
      params = "-s -d 200 -c 5,7,8 -o min,distinct,max",
      check.zero.based=F,
      check.chr=F,
      verbose=F
    );
    colnames(allbed.sort.merge) = c("chr","start","end","strand","evalue","gene","weightFold");
    allbed.sort.merge[,"region"] = paste(allbed.sort.merge[,"chr"],":",allbed.sort.merge[,"start"],"-",allbed.sort.merge[,"end"],sep="");
    result = allbed.sort.merge[,c("chr","start","end","region","evalue","strand","gene","weightFold")]
  } else {
    result = allbed.sort[,c("chr","start","end","region","evalue","strand","gene","weightFold")]
  }
}

write.table(result,file=outfile,col.names=T,row.names=F,sep="\t",quote=F);

