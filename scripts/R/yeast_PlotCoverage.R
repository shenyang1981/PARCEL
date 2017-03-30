## adding filters to the called regions
# a). p-val < 10/n
# b). rel. abundance within transcript > median+sd
# c). abs(fold change) > 2
# d). normalized counts	 of two conditions > 10
# e). both test conditions pass mean+2sd of all control	conditions

plot_candidate = function(v1all = NULL,flank = 5, myposition = NULL, id = NULL){

	tmp = v1all[(myposition[1]-flank):(myposition[2]+flank),];
	singleCol = grep(id,colnames(tmp));
	
	for(i in 1:dim(tmp)[2]){
		x = 1:dim(tmp)[1];
		y = as.numeric(tmp[,i]);
		ylim=range(tmp);
		if(i == 1){
			plot(x,y,ylim=ylim,type="l",col="grey");
		}
		if(sum(singleCol==i)>0){
			lines(x,y,ylim=ylim,col="red");
		}
	}
}



#setwd("/mnt/projects/sunm/others/yeast/analysis_combined")
source("/mnt/projects/sunm/bacteria/useful_Rscript.R")
load("/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata")


args = commandArgs(T);
args = args[-1];
print(args);

#rowinfotab = "/mnt/projects/sunm/others/yeast/fastq_combined/atp_1.tab"
#covcutoff = 10
#conditions = "tpp"
#resultdir = "analysis_combined/";
#covinfo = "yeast_combined_v1all.Rdata";
#sffile = "yeast_combined_v1all_edgeR_sf.Rdata";

rowinfotab = args[1]; # "/mnt/projects/sunm/others/yeast/fastq_combined/apt_1.tab"
covcutoff = args[2]; # 10
conditions = args[3]; # atp
resultdir = args[4]; # analysis_combined/
covinfo = args[5]; # "yeast_combined_v1all.Rdata"
sffile = args[6]; # "yeast_combined_v1all_edgeR_sf.Rdata"

v1 <- readtab(rowinfotab)
l <- sapply(v1,length)+200
p1 <- as.character(unlist(mapply(rep,names(l),l)))
p2 <- as.character(unlist(lapply(l,function(x) c(1:x))))

load(paste(resultdir,covinfo,sep=""));
n <- length(which(rowSums(v1all)>NCOL(v1all))) # average count > 1 for given position
thr <- 10/n # pvalue threshold 
v1names <- rownames(v1all)

normalizedCountFilter = 10;
numConditionPassNCF = 2;

load(paste(resultdir,sffile,sep=""));

#conditions <- c("arg","asn","asp","atp","cys","glucose","gly","his",
#                "lys","oxal","pro","sam","ser","tpp","trp","tyr","utp","val")

v1all_scaled = t(t(v1all)*sf);

result = read.table(paste(resultdir,"yeast_combined_",conditions,"_output2_wfilters.txt",sep=""),header=T,sep="\t");
pdf(paste(resultdir,"yeast_combined_",conditions,"_output2_coveragePlot.pdf",sep=""));
for(i in 1:dim(result)[1]){
	mypos = c(result[i,"begin"],result[i,"end"]);
	plot_candidate(v1all = v1all_scaled,myposition = mypos,id=conditions);
}
dev.off();
	
  
