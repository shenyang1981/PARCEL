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


rowinfotab = args[1]; # "/mnt/projects/sunm/others/yeast/fastq_combined/apt_1.tab"
covcutoff = as.numeric(args[2]); # 10
conditions = args[3]; # atp
resultdir = args[4]; # analysis_combined/
covinfo = args[5]; # "yeast_combined_v1all.Rdata"
sffile = args[6]; # "yeast_combined_v1all_edgeR_sf.Rdata"
gStruc = args[7]; # "/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata"
sampleinfofile = args[8]; # sampleinfo.txt
batchid = args[9];
sampleInfo = read.table(sampleinfofile,header=T,sep="\t",quote='"',stringsAsFactor=F);
sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == batchid,c("Condition","LibID")];

if(grepl("Rdata$",gStruc)){
	load(gStruc);
} else {
	gStruc = read.table(gStruc,header=T,sep="\t",quote='"');
}


v1 <- readtab(rowinfotab)
l <- sapply(v1,length)+200
indexbyid = cumsum(l); # get the row index for each gene
#p1 <- as.character(unlist(mapply(rep,names(l),l)))
#p2 <- as.character(unlist(lapply(l,function(x) c(1:x))))

#load(paste(resultdir,covinfo,sep=""));
load(covinfo);
n <- length(which(rowSums(v1all)>NCOL(v1all))) # average count > 1 for given position
thr <- 10/n # pvalue threshold 
v1names <- rownames(v1all)

normalizedCountFilter = 10;
numConditionPassNCF = 2;
foldcutoff = 2;
#load(paste(resultdir,sffile,sep=""));
load(sffile);

#conditions <- c("arg","asn","asp","atp","cys","glucose","gly","his",
#                "lys","oxal","pro","sam","ser","tpp","trp","tyr","utp","val")

#v1all_scaled = t(t(v1all)*sf);
#browser();
for (c in conditions){
  print(c)
  load(paste(resultdir,"fastq2_",c,"_output10.Rdata",sep="")) # output
  load(paste(resultdir,"etTable_",c,".Rdata",sep="")) # etTable
  group_names <- rep(FALSE,NCOL(v1all))
  conditionLibIDs = sampleInfo[sampleInfo[,"Condition"]==c,"LibID"];
  group_names[!is.na(match(colnames(v1all),conditionLibIDs))] <- TRUE
  
  fc <- maxpos <- dist <- min2 <- pcounts <- lowpass <- numeric()
  pcpos <- lppos <- character()
  #output = output[output[,"geneID"]=="orf19.6854",];
  if(NROW(output)>0){
    pvaluefilter = rep("nopass",NROW(output)); # adjusted pvalue <= 0.1
    covfilter = rep("nopass",NROW(output)); # the sites with significant changes should has enough site coverage within candidate region
    foldfilter =  rep("nopass",NROW(output)); # fold change >=1 1 for at least one site
    sitecovfilter = rep("nopass",NROW(output)); # coverage filter for sites with engouh fold change 
    changefilter = rep("nopass",NROW(output)); # read coverage of treatment samples should be consistently greater/less than read coverage of other control samples
    
  for (i in 1:NROW(output)){
    #cat(".");
    print(output[i,7]);
    # relative abundance for each transcript
#   	if(i != 732){
#       next;
#   	}
#     browser();
    if(output[i,"geneID"]=="orf19.6854"){
      #browser();
    } else {
      #next;
    }

    #gpos <- rownames(etTable)%in%(v1names[p1==output[i,7]])
    upboundinV1 = indexbyid[as.character(output[i,7])] - l[as.character(output[i,7])] + 1;
    lowerboundinV1 = indexbyid[as.character(output[i,7])];
    gpos <- rownames(etTable)%in%(v1names[upboundinV1:lowerboundinV1])
    et <- etTable[gpos,]  
    
    #logcpm et[,2], min2[i] is the relative abundance. median + sd
#    min2[i] <- median(et[is.finite(et[[2]]),2])+sd(et[is.finite(et[[2]]),2]) this criteria is too strigent so that orf19.6854 can not pass the fold change criteria
    #min2[i] <- median(et[is.finite(et[[2]]),2])
    min2[i] = quantile(et[is.finite(et[[2]]),2],0.65)
    pvaluefilter[i] = ifelse(sum(et[[3]]<thr)>0,"pass","nopass");
		covfilter[i] = ifelse(sum(et[[2]]>=min2[i]&et[[3]]<thr)>0,"pass","nopass");
		passed <- et[et[[2]]>=min2[i]&et[[3]]<thr,]
	  passed <- passed[rownames(passed)%in%v1names[output[i,1]:output[i,2]],]
    ## passed - <- positions passed min2, thr,

		#gStruc is the CDS information start, end of CDS, which is used to annotate the biotype of candidates.
    g <- gStruc[gStruc$gene==as.character(output[i,7]),]
    fc[i] <- maxpos[i] <- dist[i] <- pcounts[i] <- lowpass[i] <- 0  
    pcpos[i] <- lppos[i] <- "none"
    
    if (NROW(passed)>0){ # pass the gene coverage information    
      #print(i);
      fc[i] <- 2^(max(abs(passed[[1]]))) # fold change filter
      #if (fc[i]>2){ # if fold change > 2fold
      if (fc[i]>foldcutoff){ # if fold change > 2fold
					foldfilter[i] = "pass";				
#         if(length(maxpos[i]) != length(as.numeric(p2[v1names==rownames(passed[which.max(abs(passed[[1]])),])]))){
#           browser();
#         }
        #maxpos[i] <- as.numeric(p2[v1names==rownames(passed[which.max(abs(passed[[1]])),])])[1]
        maxpos[i] = as.numeric(substring(rownames(passed[which.max(abs(passed[[1]])),]),nchar(as.character(output[i,7]))+1)); # the rowname is like geneid+pos, so it can be extracted by substr()

        j <- maxpos[i]
        
        # get the relative position in transcript, 5UTR negative, CDS position, 3UTR positive + 100.
        if (j<=g[[2]]) {
          dist[i] <- j-g[[2]]-1
        } else if (j>=g[[3]]){
          dist[i] <- j-g[[3]]+1001
        } else dist[i] <- (j-g[[2]]+1)*1000/(g[[3]]-g[[2]]-1)
        passed_pos <- which(v1names%in%rownames(passed[abs(passed[[1]])>=log2(foldcutoff),])) ## only those positions that also pass FC filters
				
				#
        if (length(passed_pos)>1){
          passed_counts <- t(t(v1all[passed_pos,])*sf) # sf is from  yeast_combined_v1all_edgeR_sf.Rdata
          pcounts[i] <- as.numeric(length(which(rowSums(passed_counts>normalizedCountFilter)>=numConditionPassNCF))>0)
	        if(pcounts[i]>0) {
            pcpos[i] <- paste(passed_pos[which(rowSums(passed_counts>normalizedCountFilter)>=numConditionPassNCF)],collapse=";")
					  sitecovfilter[i] = "pass";
	        }
          passed2 <- passed[match(rownames(passed_counts),rownames(passed)),]  	     
	        sign <- ifelse(passed2[[1]]>0,TRUE,FALSE) ## TRUE for up, FALSE for down

	        lpstatus <- unlist(lapply(1:NROW(passed2),function(x){
  		      if (sign[x]){
  		        return(min(passed_counts[x,group_names])>(mean(passed_counts[x,!group_names])+2*sd(passed_counts[x,!group_names])))
  		      } else return(max(passed_counts[x,group_names])<(mean(passed_counts[x,!group_names])-2*sd(passed_counts[x,!group_names])))
	        }));

	        lowpass[i] <- sum(as.numeric(lpstatus))>0
	        if (lowpass[i]>0) lppos[i] <- paste(passed_pos[lpstatus],collapse=";")

        } else {
          passed_counts <- v1all[passed_pos,]*sf
          pcounts[i] <- as.numeric(sum(passed_counts>normalizedCountFilter)>=numConditionPassNCF)
	        if (pcounts[i]>0) {
            pcpos[i] <- passed_pos;
            sitecovfilter[i] = "pass";
	        }

          if (passed$logFC[abs(passed[[1]])>log2(foldcutoff)]>0) { # increasing/decreasing of V1 signal 
            lowpass[i] <- as.numeric(min(passed_counts[group_names]) > mean(passed_counts[!group_names])+2*sd(passed_counts[!group_names]))
          } else {
            lowpass[i] <- as.numeric(max(passed_counts[group_names]) < mean(passed_counts[!group_names])-2*sd(passed_counts[!group_names]))            
          }
	        if (lowpass[i]>0) {
	        	lppos[i] <- passed_pos;
	        	changefilter[i] = "pass";
	        }
        }
      }
    }
  }
 }
  originaloutput <- data.frame(output,maxPos=maxpos,maxFC=fc,m1sd=min2,dist=dist,pcounts,pcpos,lowpass,lppos,pvaluefilter,covfilter,foldfilter,sitecovfilter,changefilter);
  output2 <- data.frame(output,maxPos=maxpos,maxFC=fc,m1sd=min2,dist=dist,pcounts,pcpos,lowpass,lppos)[fc>foldcutoff,]
  print(paste(c,": ",NROW(output2[output2$pcounts>0&output2$lowpass>0,]),"positive",sep=""));
  if (NROW(output2[output2$pcounts>0&output2$lowpass>0,])>0){ 
     output2b <- output2[output2$pcounts>0&output2$lowpass>0,] 
     save(output2b,originaloutput,file=paste(resultdir,"combined_",c,"_output2_wfilters.Rdata",sep=""))
     write.table(output2b,file=paste(resultdir,"combined_",c,"_output2_wfilters.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
  } else {
     output2b = output2
     save(output2,originaloutput,file=paste(resultdir,"combined_",c,"_output2_wfilters.Rdata",sep=""))
     write.table(output2b,file=paste(resultdir,"combined_",c,"_output2_wfilters.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
  }
}

  
