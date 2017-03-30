args = commandArgs(T);

args = args[-1];

# args[1] cutoff
# args[2] condition
# args[3] scaling factor, 
outdir = args[1]; # ~/prog/parcel/analysis_combined/"
cutoff = as.numeric(args[2]);
sampleid = args[3];
countinfo = args[4];
sffile = paste(outdir,args[5],sep=""); # yeast_combined_v1all_edgeR_sf.Rdata

outdir = "analysis_combined_yeastall/"; 
cutoff = as.numeric(10);
sampleid = "pro";
countinfo = "yeast_combined_v1all.Rdata";
sffile = paste(outdir,"yeast_combined_v1all_edgeR_sf.Rdata",sep=""); # yeast_combined_v1all_edgeR_sf.Rdata



#inputdir = "~/prog/parcel/analysis_combined/";

load(paste(outdir,countinfo,sep=""));
require(edgeR)
countable <- v1all[rowSums(v1all)>cutoff,]

i = sampleid;
group_names <- rep("controls",NCOL(countable))
group_names[grep(i,colnames(countable))] <- i

colnames(countable) = group_names;
#rownames(countdata) = as.character(countdataAll[,"ID"]);


y <- DGEList(counts = countable,group=group_names,genes=data.frame(Length=rep(1,NROW(countable))))
y <- calcNormFactors(y)
if(!file.exists(sffile)){
	sf <- y$samples$norm.factors
	save(sf,file=sffile)
}
y <- estimateCommonDisp(y, verbose = TRUE)
y <- estimateTagwiseDisp(y)
et <- exactTest(y, dispersion = "auto")
etTable <- et$table
save(etTable,file=paste(outdir,"combined_v1all_etTable_",i,".Rdata",sep=""))






coldata = data.frame("shRNA" = c("shc","sh10","shc","sh10"),"chip"=c("aff3","aff3","input","input"));

rownames(coldata) = colnames(countdata);

dds = DESeqDataSetFromMatrix(countData = countdata,colData = coldata,design = ~ shRNA + chip);

dds = dds[rowSums(counts(dds))>1,];

dds = estimateSizeFactors(dds);
dds <- estimateDispersions(dds)
# rld = rlog(dds,blind=F);
# plot(log2(counts(dds,normalized=T)[,1:2]+1),pch=16,cex=0.3);
# plot(assay(rld)[,1:2],pch=16,cex=0.3);
# abline(0,1)
# abline(1,1)
# abline(-1,1);

#plotPCA(rld,intgroup = c("shRNA","chip"))

dds = DESeq(dds,test="LRT",reduced=~ chip);
result = results(dds);

