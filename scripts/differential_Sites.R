suppressPackageStartupMessages(library("data.table"))
args = commandArgs(T);

args = args[-1];

# args[1] cutoff
# args[2] condition
# args[3] scaling factor, 
outdir = args[1]; # ~/prog/parcel/analysis_combined/"
cutoff = as.numeric(args[2]);
sampleid = args[3];
controlid = args[4];
vsall = as.logical(args[5]);
sampleinfofile = args[6];
countinfo = args[7];
sffile = args[8]; # yeast_combined_v1all_edgeR_sf.Rdata
batchid = args[9]; # comparison batch
ismerge = ifelse(!is.null(args[10]) & args[10]=="T",T,F);
#browser();
sampleInfo = read.table(sampleinfofile,header=T,sep="\t",quote='"',stringsAsFactor=F);

sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == batchid,];

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
#browser();
#inputdir = "~/prog/parcel/analysis_combined/";

conditions = sampleInfo[,"Condition"];

#v1all = v1all[,sampleInfo[,"LibID"],with=F];

suppressPackageStartupMessages(require(edgeR));
countable <- as.matrix(v1all[rowSums(v1all[,-c("chr","pos"),with=F])>=cutoff,-c("chr","pos"),with=F]);
#colnames(countable) = conditions;

#conditions <- c("ala","arg","asn","asp","atp","b12","camp","ctp","cys","fmn")
# conditions <- c("g6p","gln","glu","glucose","gly","gtp","his","iso","isoa","leu")
# conditions <- c("lys","mala","met","nadh","nasu","oxal","phe","pro","pyr","sam")
# conditions <- c("ser","sper","thr","tpp","trp","tyr","utp","val")
i = sampleid;
group_names <- rep("controls",NCOL(countable))
#group_names[grepl(i,conditions)] <- i
group_names[!is.na(match(conditions,i))] <- i
groupLevel = factor(group_names,levels = c("controls",i));

#y <- DGEList(counts = countable,group=groupLevel,genes=data.frame(Length=rep(1,NROW(countable))))
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
save(etTable,file=paste(outdir,"etTable_",i,".Rdata",sep=""))
save(v1all,file=paste(outdir,"covinfo_",i,".Rdata",sep=""))

