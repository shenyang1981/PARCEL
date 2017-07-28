#!/usr/bin/env Rscript
setwd("~/gseq/prog/parcel/GutBacteria/parcelResultSunMiao/batch1/")
load("fastq2_metpool_output10.Rdata")
load("combined_metpool_output2_wfilters.Rdata")

output[,"chr"] = gsub("(.*?):.*","\\1",output[,"bin"])
output = output[output[,"pstart"]!=output[,"pend"],]
output.sort = bedr.sort.region(output[,c("chr","pstart","pend","bin")],check.zero.based = F,check.chr = F,check.valid = F,check.merge = F,engine = "bedtools")
output.merge = bedr(
  engine = "bedtools",
  method="merge",
  input = list(i=output.sort),
  params = "-d 200 -c 4 -o collapse",check.chr = F,check.zero.based = F,verbose = F,check.merge=F,check.valid=F,check.sort=F);
output.merge[,"newbin"] = paste(output.merge[,1],":",output.merge[,2],"-",output.merge[,3],sep="");
rownames(output.merge) = NULL;
mergebins = unlist(apply(output.merge[,c("newbin","bin")],1,function(x){y = unlist(strsplit(x[2],split=',')); result=rep(x[1],length(y));names(result) = y;return(result)}))
output[,"newbin"] = mergebins[output[,"bin"]];

ribobed = read.table("~/gseq/prog/parcel/RegPrecise/allRNA/riboMergedHMPREFGut.sort.bed",header=F,sep="\t")
ribobed.ext = ribobed
ribobed.ext[,3] = ribobed.ext[,3]+200
ribobed.ext[,2] = ifelse(ribobed.ext[,2]-200 >= 1, ribobed.ext[,2]-200 ,1)

colnames(output.sort)[1:3] = c("chr","start","end")
colnames(ribobed.ext)[1:3] = c("chr","start","end")
outputInRibo = bedInterSect(region1 = output.sort,region2=ribobed.ext)

# originaloutputRibo = merge(originaloutput,outputInRibo[,c("bin","V4")],by="bin",all.x=T)
# colnames(originaloutputRibo)[dim(originaloutputRibo)[2]] = "KnownRibo";

allcols = colnames(originaloutput);
filtercols = allcols[grepl("filter",allcols)];

outputRibo = merge(output[,c("bin","newbin","Evalue")],outputInRibo[,c("bin","V4")],by="bin",all.x=T)
colnames(outputRibo)[dim(outputRibo)[2]] = "KnownRibo";
outputRibo = merge(outputRibo,originaloutput[,c("bin",filtercols)],by="bin",all.x=T);
outputRibo[outputRibo[,"Evalue"]<=5,"Evaluefilter"] = "pass";


allcols = colnames(outputRibo);
filtercols = allcols[grepl("filter",allcols)];
filters = c("Evaluefilter","pvaluefilter","covfilter","foldfilter","sitecovfilter","changefilter");

knownriboBin = tapply(outputRibo[,"KnownRibo"],mergebins[outputRibo[,"bin"]],function(x){sum(!is.na(x))>0;});
newbinResult = apply(outputRibo[,c(filters)],2,function(x){
  result = tapply(x,mergebins[outputRibo[,"bin"]],function(x){sum(!is.na(x) & x=="pass")>0;})
  return(unlist(result));});
newbinResult = cbind(knownriboBin,newbinResult);



#all = apply(outputRibo[,filtercols],2,function(x){table(!is.na(x) & x=="pass",!is.na(outputRibo[,"KnownRibo"]))});
all = apply(newbinResult[,-1],2,function(x){table(x,newbinResult[,1])});
all = all[c(2,4),];
all = cbind(table(!is.na(outputRibo[,"KnownRibo"])),all);
colnames(all)[1] = "All";


xlabs = gsub("(.*)filter","\\1",filters);
barplot(all[1,filters],names.arg = xlabs,col=trop[2])
#originaloutputRibo[,"filterstr"] = apply(originaloutputRibo[,filtercols],1,function(x){paste(ifelse(x=="pass","P","NP"),collapse="_")});
#tmp = as.matrix(table(is.na(originaloutputRibo[,"KnownRibo"]),originaloutputRibo[,"filterstr"]));
#tmp[,1] = rowSums(tmp);
#colnames(tmp) = c("All",filtercols[-1])