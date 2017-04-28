
# convert tab coverage into data.frame v1all
# row is chromsome:position and column is library id.

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
dir.create(outdir);

tabdir = args[2]; #/mnt/projects/sunm/others/yeast/fastq_combined/
candsfile = args[3]; #"~/prog/parcel/yeast_candidates.txt"
countinfo = args[4]; # yeast_combined_v1all.Rdata

v1all <- NULL

# load candidates

cands = read.table(candsfile,header=F,sep="\t",quote='"');

for (condition in as.character(cands[,1])){
  v1 <- readtab(paste(tabdir,"/",condition,".tab",sep=""))
  v1 <- lapply(v1,function(x) c(x,rep(0,200)))
  yeast <- do.call("c",v1)
  v1all <- cbind(v1all,yeast)      
}
colnames(v1all) <- as.character(cands[,1])         

#save(v1all,file=paste(outdir,countinfo,sep=""));
save(v1all,file=countinfo);


