
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

sampleInfo = read.table(sampleinfofile,header=T,sep="\t",quote='"',stringsAsFactor=F);

sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == batchid,c("Condition","LibID")];

if(!vsall){
  sampleInfo = sampleInfo[!is.na(match(sampleInfo[,"Condition"],c(sampleid,controlid))),];
}
#browser();
#inputdir = "~/prog/parcel/analysis_combined/";

# load v1all from file countinfo
load(countinfo);

sampleInfo = sampleInfo[!is.na(match(sampleInfo[,"LibID"],colnames(v1all))),];

conditions = sampleInfo[,"Condition"];
v1all = v1all[,sampleInfo[,"LibID"]];

require(edgeR)
countable <- v1all[rowSums(v1all)>cutoff,]
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

y <- DGEList(counts = countable,group=groupLevel,genes=data.frame(Length=rep(1,NROW(countable))))
y <- calcNormFactors(y)
if(!file.exists(sffile)){
  sf <- y$samples$norm.factors
	save(sf,file=sffile)
}
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y, dispersion = "auto")
etTable <- et$table
save(etTable,file=paste(outdir,"etTable_",i,".Rdata",sep=""))
save(v1all,file=paste(outdir,"covinfo_",i,".Rdata",sep=""))

