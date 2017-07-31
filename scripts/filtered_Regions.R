#!/usr/bin/env Rscript
## adding filters to the called regions
# a). adjusted p-val < 0.1
# b). rel. abundance within region +- 120bp > median of abundance of region +- 120bp
# c). abs(fold change) > 2
# d). normalized counts	 of two conditions > 10
# e). both test conditions pass mean+2sd of all control	conditions

suppressPackageStartupMessages(library("data.table"));

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

filterOutput = function(x,thr=NULL,group_names=NULL){
  print(x[["geneID"]]);
  
  # get the fold change and pvalue
  et = etTable[chr==x[["chr"]],]
  
  # logcpm et[,2], min2[i] is the relative abundance. median + sd
  localcovfilter = 0.5;
  extendLen = 120; # how long is the length for extension
  tmpstart = as.numeric(as.character(x[["pstart"]])) - extendLen;
  tmpend = as.numeric(as.character(x[["pend"]])) + extendLen;
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
  
  return(
    list( maxPos= as.numeric(maxpos), maxFC=fc, m1sd=min2, dist=dist, pcounts = pcounts,
          pcpos = pcpos,lowpass = lowpass,lppos = lppos,pvaluefilter = pvaluefilter,
          covfilter = covfilter,foldfilter = foldfilter,sitecovfilter = sitecovfilter, changefilter = changefilter));
    
}


args = commandArgs(T);
args = args[-1];
print(args);

##########################################################################
# Note, to run this pipeline, at least 2 control samples are required
##########################################################################


covcutoff = as.numeric(args[1]); # coverage cutoff
conditions = args[2]; # treatment
resultdir = args[3]; # 
covinfo = args[4]; # coverage information from "combined_v1all.Rdata"
sffile = args[5]; # scaling factor to calculate normalized counts,  "combined_v1all_edgeR_sf.Rdata"
gStruc = args[6]; # CDS information 
sampleinfofile = args[7]; # sampleinfo.txt
batchid = args[8]; # batchID
ismerge =  ifelse(!is.null(args[9]) & args[9]=="T",T,F); # whether merge libraries from same condition

# process sample information ---------------------------
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
n = sum(rowSums(v1all[,c(-1,-2),with=F]) > (NCOL(v1all)-2)); 
thr = 0.1; # pvalue threshold 


normalizedCountFilter = 10;
numConditionPassNCF = 2;
foldcutoff = 2;


usedCol = sampleInfo[,"LibID"];  
group_names <- rep(FALSE,NROW(sampleInfo))
names(group_names) = usedCol;
conditionLibIDs = sampleInfo[sampleInfo[,"Condition"]==conditions,"LibID"];
group_names[!is.na(match(names(group_names),conditionLibIDs))] <- TRUE

  
if(NROW(output)>0){
  mainCols = names(output);
  output[,evaluefilter:="nopass"];
  output[E_value<=10,evaluefilter:="pass"];
  v1all  = v1all[chr %in% output[evaluefilter=="pass",chr],];
  output[,rowid:=.I];
  returnVars = c("maxPos","maxFC","m1sd","dist","pcounts","pcpos","lowpass","lppos","pvaluefilter","covfilter","foldfilter","sitecovfilter","changefilter")
  output[evaluefilter == "pass", (returnVars):= filterOutput(.SD,thr=thr,group_names=group_names),by = rowid];
  originaloutput <- output;
  finalCols = c(mainCols,"maxPos","maxFC","m1sd","dist","pcounts","pcpos","lowpass","lppos");
  output2 <- output[changefilter=="pass",finalCols,with=F];
  print(paste(conditions,": ",NROW(output2),"positive",sep=""));
  save(output2,originaloutput,file=paste(resultdir,"combined_",conditions,"_output2_wfilters.Rdata",sep=""))
  write.table(output2,file=paste(resultdir,"combined_",conditions,"_output2_wfilters.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
}
