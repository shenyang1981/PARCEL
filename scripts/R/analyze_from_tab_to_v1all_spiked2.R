readtab <- function(filename){
  v1 <- scan(filename,what="",sep="\n")
  v1 <- strsplit(v1, "\t")
  names(v1) <- sapply(v1, `[[`, 1)
  v1 <- lapply(v1, `[`, -1)
  v1 <- lapply(v1,function(x) as.numeric(strsplit(x,";")[[1]]))
  return(v1)}

args = commandArgs(T);

args = args[-1];


outdir = args[1]; #"~/prog/parcel/analysis_combined/";
tabdir = args[2]; #/mnt/projects/sunm/others/yeast/fastq_combined/
candsfile = args[3]; #"~/prog/parcel/yeast_candidates.txt"
countinfo = args[4]; # yeast_combined_v1all.Rdata

#setwd("")
v1all <- NULL
out <- read.table(paste(tabdir,"mapping_combined.out",sep=""));
# load candidates

sampleid = gsub("(.*)_(.*)","\\1",out[,1],perl=T);
cands = read.table(candsfile,header=F,sep="\t",quote='"');
out = out[!is.na(match(sampleid,cands[,1])),];


for (condition in as.character(out[,1])){
  v1 <- readtab(paste(tabdir,condition,".tab",sep=""))
  v1 <- lapply(v1,function(x) c(x,rep(0,200)))
  yeast <- do.call("c",v1)
  v1all <- cbind(v1all,yeast)      
}
colnames(v1all) <- as.character(out[,1])         

save(v1all,file=paste(outdir,countinfo,sep=""));


# vcor <- cor(v1all)
# save(vcor,file="yeast_combined_v1all_cor.Rdata")
# write.table(vcor,file="yeast_combined_v1all_cor.txt",sep="\t",quote=FALSE)

# vcor <- cor(log(v1all+1))
# save(vcor,file="yeast_combined_v1all_cor_log+1.Rdata")
# write.table(vcor,file="yeast_combined_v1all_cor_log+1.txt",sep="\t",quote=FALSE)

# library(corrplot)
# corrplot(vcor,method="circle")
# atp_2, glucose_1, h2o_3, h2o_4, perhaps sam_2

