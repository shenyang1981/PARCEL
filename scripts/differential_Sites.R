#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(edgeR));
suppressPackageStartupMessages(library("data.table"))
args = commandArgs(T);
args = args[-1];

outdir = args[1]; # ~/prog/parcel/analysis_combined/"
cutoff = as.numeric(args[2]); # coverage cutoff
sampleid = args[3]; #name of treatment
controlid = args[4]; #name of control 
vsall = as.logical(args[5]); # treatments vs. control or treatments vs. others
sampleinfofile = args[6]; # sample information
countinfo = args[7]; # coverage information stored in Rdata
sffile = args[8]; # scaling factors stored in combined_v1all_edgeR_sf.Rdata
batchid = args[9]; # batchID defined in sample information file

ismerge = ifelse(!is.null(args[10]) & args[10]=="T",T,F);
sampleInfo = read.table(sampleinfofile,header=T,sep="\t",quote='"',stringsAsFactor=F);
sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == batchid,];

# merge libraries with same condition
if(ismerge){
  sampleInfo[,"LibID"] = paste(sampleInfo[,"Condition"],sampleInfo[,"Replicates"],sep="__");
}
sampleInfo = unique(sampleInfo[,c("Condition","LibID")]);

# load v1all from file countinfo
load(countinfo);

if(!vsall){
  	sampleInfo = sampleInfo[!is.na(match(sampleInfo[,"Condition"],c(sampleid,controlid))),];
	v1all = v1all[,c("chr","pos",sampleInfo[,"LibID"]),with=F];
}

conditions = sampleInfo[,"Condition"];

# Detect differential sites using edgeR ----------------------------

countable <- as.matrix(v1all[rowSums(v1all[,-c("chr","pos"),with=F])>=cutoff,-c("chr","pos"),with=F]);

group_names <- rep("controls",NCOL(countable))
group_names[!is.na(match(conditions,sampleid))] <- sampleid
groupLevel = factor(group_names,levels = c("controls",sampleid));

y <- DGEList(counts = countable,group=groupLevel)
y <- calcNormFactors(y)
if(!file.exists(sffile)){
  sf <- y$samples$norm.factors;
  names(sf) = rownames(y$samples);
	save(sf,file=sffile)
}
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y, dispersion = "auto")
etTable <- data.table(v1all[rowSums(v1all[,-c("chr","pos"),with=F])>=cutoff,c("chr","pos"),with=F],et$table);
save(etTable,file=paste(outdir,"etTable_",sampleid,".Rdata",sep=""))
save(v1all,file=paste(outdir,"covinfo_",sampleid,".Rdata",sep=""))

