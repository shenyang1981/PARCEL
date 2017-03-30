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

#inputdir = "~/prog/parcel/analysis_combined/";

load(paste(outdir,countinfo,sep=""));
require(edgeR)
countable <- v1all[rowSums(v1all)>cutoff,]

conditions <- c("ala","arg","asn","asp","atp","b12","camp","ctp","cys","fmn")
# conditions <- c("g6p","gln","glu","glucose","gly","gtp","his","iso","isoa","leu")
# conditions <- c("lys","mala","met","nadh","nasu","oxal","phe","pro","pyr","sam")
# conditions <- c("ser","sper","thr","tpp","trp","tyr","utp","val")
#for (i in conditions){
		i = sampleid;
    group_names <- rep("controls",NCOL(countable))
    group_names[grep(i,colnames(countable))] <- i
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
#}
