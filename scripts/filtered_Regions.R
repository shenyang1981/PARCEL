## adding filters to the called regions
# a). p-val < 10/n
# b). rel. abundance within transcript > median+sd
# c). abs(fold change) > 2
# d). normalized counts	 of two conditions > 10
# e). both test conditions pass mean+2sd of all control	conditions

plot_candidate = function(v1all = NULL,flank = 5, myposition = NULL, id = NULL){

	tmp = v1all[myposition[1]:myposition[2],];
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

filterOutput = function(x){
  print(x[["geneID"]]);
  
  # get the fold change and pvalue
  #et = etTable[chr==x[["chr"]] & pos >= (x[["start"]]-200) & pos <= (x[["end"]]+200),]
  et = etTable[chr==x[["chr"]],]
  
  # logcpm et[,2], min2[i] is the relative abundance. median + sd
  # min2[i] <- median(et[is.finite(et[[2]]),2])+sd(et[is.finite(et[[2]]),2]) this criteria is too strigent so that orf19.6854 can not pass the fold change criteria
  #min2[i] <- median(et[is.finite(et[[2]]),2])
  localcovfilter = 0.5;
  extendLen = 120; # how long is the length for extension
  tmpstart = as.numeric(as.character(x[["pstart"]])) - extendLen;
  tmpend = as.numeric(as.character(x[["pend"]])) + extendLen;
  #min2 = et[is.finite(logCPM),quantile(logCPM,localcovfilter)]; # coverage cutoff for given gene
  #min2 = et[is.finite(logCPM) & pos >= tmpstart & pos <= tmpend,quantile(logCPM,localcovfilter)]; # coverage cutoff for given gene
  min2 = et[is.finite(logCPM) & pos >= tmpstart & pos <= tmpend,mean(logCPM)]; # coverage cutoff for given gene
  pvaluefilter = ifelse(sum(et[ pos >= as.numeric(as.character(x[["pstart"]])) & pos <= as.numeric(as.character(x[["pend"]])),PAdjust]<thr)>0,"pass","nopass");
  covfilter = ifelse(
                sum(
                  et[,pos] >= as.numeric(as.character(x[["pstart"]])) & et[,pos] <= as.numeric(as.character(x[["pend"]])) 
                  & et[,logCPM]>=min2 & et[,PAdjust]<thr)>0,"pass","nopass");
  passed <- et[logCPM>=min2 & PAdjust<thr & pos >= as.numeric(as.character(x[["pstart"]])) & pos <= as.numeric(as.character(x[["pend"]])),]
  foldfilter = "nopass";
  sitecovfilter = "nopass";
  changefilter = "nopass";
  ## passed - <- positions passed min2, thr,
  
  #gStruc is the CDS information start, end of CDS, which is used to annotate the biotype of candidates.
  g <- gStruc[gStruc$gene==gsub("(.*):(.*)","\\1",x[["chr"]]),]
  fc <- maxpos <- dist <- pcounts <- lowpass <- as.numeric(0);
  pcpos <- lppos <- "none"
  
  if (NROW(passed)>0){ # pass the gene coverage information    
    #print(i);
    fc <- 2^(max(abs(passed[,logFC]))) # fold change filter
    #if (fc[i]>2){ # if fold change > 2fold
    if (fc>foldcutoff){ # if fold change > 2fold
      foldfilter = "pass";				
      maxpos = passed[which.max(abs(logFC)),pos]; # get the position of site with maximum fold change
      # get the relative position in transcript, 5UTR negative, CDS position, 3UTR positive + 100.
      if (maxpos<=g[[2]]) {
        dist <- maxpos-g[[2]]-1 # 5'UTR
      } else if (maxpos>=g[[3]]){
        dist <- maxpos-g[[3]]+1001 # 3'UTR
      } else {
        dist <- (maxpos-g[[2]]+1)*1000/(g[[3]]-g[[2]]-1); # CDS scaled as 1000bp
      }
      passed_pos <- passed[abs(logFC) >= log2(foldcutoff),pos] ## only those positions that also pass FC filters
      
      if (length(passed_pos)>1){ # multiple sites passed
        passed_counts = v1all[chr==x[["chr"]] & pos %in% passed_pos,];
        norcounts = t(t(passed_counts[,usedCol,with=F])*sf[usedCol]); # sf is from  yeast_combined_v1all_edgeR_sf.Rdata
        passedIndex = apply(norcounts,1,function(x){sum( x > normalizedCountFilter) >= numConditionPassNCF});
        
        if(sum(passedIndex)>0) {
          pcpos <- paste(passed_counts[passedIndex,pos],collapse=";");
          sitecovfilter = "pass";
        }
        passed2 <- passed[pos %in% passed_counts[,pos],]  	     
        sign <- ifelse(passed2[,logFC]>0,TRUE,FALSE) ## TRUE for up, FALSE for down
        
        lpstatus <- unlist(
          lapply(1:NROW(passed2),function(x){
            if (sign[x]){
              return(min(norcounts[x,group_names])>(mean(norcounts[x,!group_names])+2*sd(norcounts[x,!group_names])))
            } else {
              return(max(norcounts[x,group_names])<(mean(norcounts[x,!group_names])-2*sd(norcounts[x,!group_names])))
            }
          }));
        lowpass <- as.numeric(sum(as.numeric(lpstatus))>0);
        if (lowpass>0) {
          lppos <- paste(passed_pos[lpstatus],collapse=";")
          changefilter = "pass";
        }
      } else {
        passed_counts <- v1all[chr==x[["chr"]] & pos == passed_pos,];
        norcounts = unlist(passed_counts[,usedCol,with=F])*sf[usedCol];
        pcounts <- as.numeric(sum(norcounts>normalizedCountFilter)>=numConditionPassNCF)
        if (pcounts>0) {
          pcpos <- as.character(passed_pos);
          sitecovfilter = "pass";
        }
        if (passed[pos == passed_pos,logFC]>0){ # increasing/decreasing of V1 signal 
          lowpass <- as.numeric(min(norcounts[group_names]) > mean(norcounts[!group_names])+2*sd(norcounts[!group_names]))
        } else {
          lowpass <- as.numeric(max(norcounts[group_names]) < mean(norcounts[!group_names])-2*sd(norcounts[!group_names]))            
        }
        if (lowpass>0) {
          lppos <- as.character(passed_pos);
          changefilter = "pass";
        }
      }
    }
  }
  
  # to use do.call(rbind.data.frame,function), the result should be list of list
  #resultlist = list();
  #resultlist[[1]] = list(maxPos=as.numeric(maxpos),maxFC=as.numeric(fc));
  
  return(
    #list(maxPos=as.numeric(maxpos),maxFC=as.numeric(fc))
    list( maxPos= as.numeric(maxpos), maxFC=fc, m1sd=min2, dist=dist, pcounts = pcounts,
          pcpos = pcpos,lowpass = lowpass,lppos = lppos,pvaluefilter = pvaluefilter,
          covfilter = covfilter,foldfilter = foldfilter,sitecovfilter = sitecovfilter, changefilter = changefilter));
    
}


#### for debug #####
# rowinfotab = "~/gseq/prog/parcel/candida_tab/h2o_1.tab"
# covcutoff = 10
# conditions = "met"
# resultdir = "./"
# covinfo = "candida_combined_v1all.Rdata"
# sffile = "edgeR_sf.Rdata"
# gStruc = "~/gseq/prog/parcel/Candida_Snyder_gStruc.txt"
# sampleinfofile = "../sampleList_Calb.txt"
# batchid = "batch1"
#######################

#setwd("/mnt/projects/sunm/others/yeast/analysis_combined")
#source("/mnt/projects/sunm/bacteria/useful_Rscript.R")
#setwd("~/gseq/prog/parcel/");
suppressPackageStartupMessages(library("data.table"));
source("~/gseq/prog/parcel/C.albicans/sunScripts/useful_Rscript.R")

args = commandArgs(T);
args = args[-1];
print(args);
##########################################################################
# Note, to run this pipeline, at least 2 control samples are required
##########################################################################


# rowinfotab = "~/gseq/prog/parcel/C.albicans/tabs/batch1/V1_2.tab"
# covcutoff = 1
# conditions = "met"
# resultdir = "~/gseq/prog/parcel/C.albicans/parcelResultSunMiao/batch1/";
# covinfo = "~/gseq/prog/parcel/C.albicans/parcelResultSunMiao/batch1/covinfo_met.Rdata";
# sffile = "~/gseq/prog/parcel/C.albicans/parcelResultSunMiao/batch1/edgeR_met_sf.Rdata";
# gStruc = "~/gseq/prog/parcel/C.albicans/sunScripts/Candida_Snyder_gStruc.txt"; # "/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata"
# sampleinfofile = "~/gseq/prog/parcel/C.albicans/sampleList_Calb.txt"
# batchid = "batch1"

# rowinfotab = "~/gseq/prog/parcel/Fly/tabs/batch100/MERGE02.tab"
# covcutoff = 1
# conditions = "glucose"
# resultdir = "~/gseq/prog/parcel/Fly/parcelResultSunMiao/batch100/";
# covinfo = "~/gseq/prog/parcel/Fly/parcelResultSunMiao/batch100/covinfo_glucose.Rdata";
# sffile = "~/gseq/prog/parcel/Fly/parcelResultSunMiao/batch100/edgeR_glucose_sf.Rdata";
# gStruc = "~/gseq/prog/database/Genome/Flybase/Dmel/FB2016_02/cdsinfo.txt"; # "/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata"
# sampleinfofile = "~/gseq/prog/parcel/Fly/sampleList_Fly.txt"
# batchid = "batch100"

covcutoff = 10
conditions = "tpp"
resultdir = "~/gseq/prog/parcel/Fly/parcelResultSunMiao/batch100/";
covinfo = "~/gseq/prog/parcel/Fly/parcelResultSunMiao/batch100/covinfo_tpp.Rdata";
sffile = "~/gseq/prog/parcel/Fly/parcelResultSunMiao/batch100/edgeR_tpp_sf.Rdata";
gStruc = "~/gseq/prog/database/Genome/Flybase/Dmel/FB2016_02/cdsinfo.txt"; # "/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata"
sampleinfofile = "~/gseq/prog/parcel/Fly/sampleList_Fly.txt"
batchid = "batch100"

covcutoff = 1
conditions = "met"
resultdir = "~/gseq/prog/parcel/C.albicans/parcelResultSunMiaoFast/batch1/";
covinfo = "~/gseq/prog/parcel/C.albicans/parcelResultSunMiaoFast/batch1/covinfo_met.Rdata";
sffile = "~/gseq/prog/parcel/C.albicans/parcelResultSunMiaoFast/batch1/edgeR_tpp_sf.Rdata";
gStruc = "~/gseq/prog/parcel/C.albicans/sunScripts/Candida_Snyder_gStruc.txt"; # "/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata"
sampleinfofile = "~/gseq/prog/parcel/C.albicans//sampleList_Calb.txt"
batchid = "batch1"
ismerge=F

covcutoff = as.numeric(args[1]); # 10
conditions = args[2]; # atp
resultdir = args[3]; # analysis_combined/
covinfo = args[4]; # "yeast_combined_v1all.Rdata"
sffile = args[5]; # "yeast_combined_v1all_edgeR_sf.Rdata"
gStruc = args[6]; # "/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata"
sampleinfofile = args[7]; # sampleinfo.txt
batchid = args[8];
ismerge =  ifelse(!is.null(args[9]) & args[9]=="T",T,F);
sampleInfo = read.table(sampleinfofile,header=T,sep="\t",quote='"',stringsAsFactor=F);
sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == batchid,];

if(ismerge){
  sampleInfo[,"LibID"] = paste(sampleInfo[,"Condition"],sampleInfo[,"Replicates"],sep="__");
}
sampleInfo = unique(sampleInfo[,c("Condition","LibID")]);

load(covinfo); # load coverage information, v1all

sampleInfo = sampleInfo[!is.na(match(sampleInfo[,"LibID"],names(v1all)[c(-1,-2)])),];

if(grepl("Rdata$",gStruc)){
	load(gStruc);
} else {
	gStruc = read.table(gStruc,header=T,sep="\t",quote='"');
}


load(sffile); # load factors for normalization, sf
load(paste(resultdir,"etTable_",conditions,".Rdata",sep="")) # load edgeR result, etTable
etTable[,PAdjust:=p.adjust(PValue,method="fdr")]
load(paste(resultdir,"fastq2_",conditions,"_output10.Rdata",sep="")) # load Evalue result, output
#thr <- 10/dim(etTable)[1]; # pvalue threshold
n = sum(rowSums(v1all[,c(-1,-2),with=F]) > (NCOL(v1all)-2)); 
#thr = 10/n; # pvalue threshold
thr = 0.1; # pvalue threshold 


normalizedCountFilter = 10;
numConditionPassNCF = 2;
foldcutoff = 2;





#browser();

usedCol = sampleInfo[,"LibID"];  
group_names <- rep(FALSE,NROW(sampleInfo))
names(group_names) = usedCol;
conditionLibIDs = sampleInfo[sampleInfo[,"Condition"]==conditions,"LibID"];
group_names[!is.na(match(names(group_names),conditionLibIDs))] <- TRUE

  
fc <- maxpos <- dist <- min2 <- pcounts <- lowpass <- numeric()
pcpos <- lppos <- character()

if(NROW(output)>0){
  mainCols = names(output);
  output[,evaluefilter:="nopass"];
  output[E_value<=10,evaluefilter:="pass"];
  v1all  = v1all[chr %in% output[evaluefilter=="pass",chr],];
  output[,rowid:=.I];
  returnVars = c("maxPos","maxFC","m1sd","dist","pcounts","pcpos","lowpass","lppos","pvaluefilter","covfilter","foldfilter","sitecovfilter","changefilter")
  output[evaluefilter == "pass", (returnVars):= filterOutput(.SD),by = rowid];
  originaloutput <- output;
  finalCols = c(mainCols,"maxPos","maxFC","m1sd","dist","pcounts","pcpos","lowpass","lppos");
  output2 <- output[changefilter=="pass",finalCols,with=F];
  print(paste(conditions,": ",NROW(output2),"positive",sep=""));
  save(output2,originaloutput,file=paste(resultdir,"combined_",conditions,"_output2_wfilters.Rdata",sep=""))
  write.table(output2,file=paste(resultdir,"combined_",conditions,"_output2_wfilters.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
}

  
