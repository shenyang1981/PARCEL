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



#setwd("/mnt/projects/sunm/others/yeast/analysis_combined")
#source("/mnt/projects/sunm/bacteria/useful_Rscript.R")
setwd("~/gseq/prog/parcel/");
source("useful_Rscript.R")

args = commandArgs(T);
args = args[-1];
print(args);

# rowinfotab = "candida_tab/h2o_1.tab"
# covcutoff = 100
# conditions = "TPP1"
# resultdir = "analysis_combined_candida/";
# covinfo = "candida_combined_v1all.Rdata";
# sffile = "candida_combined_v1all_edgeR_sf.Rdata";
# gStruc = "Candida_Snyder_gStruc.txt"; # "/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata"

rowinfotab = args[1]; # "/mnt/projects/sunm/others/yeast/fastq_combined/apt_1.tab"
covcutoff = as.numeric(args[2]); # 10
conditions = args[3]; # atp
resultdir = args[4]; # analysis_combined/
covinfo = args[5]; # "yeast_combined_v1all.Rdata"
sffile = args[6]; # "yeast_combined_v1all_edgeR_sf.Rdata"
gStruc = args[7]; # "/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata"

if(grepl("Rdata$",gStruc)){
	load(gStruc);
} else {
	gStruc = read.table(gStruc,header=T,sep="\t",quote='"');
}


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

for (c in conditions){
  print(c)
  load(paste(resultdir,"fastq2_",c,"_output10.Rdata",sep="")) # output
  load(paste(resultdir,"combined_v1all_etTable_",c,".Rdata",sep="")) # etTable
  group_names <- rep(FALSE,NCOL(v1all))
  group_names[grep(c,colnames(v1all))] <- TRUE
  
  fc <- maxpos <- dist <- min2 <- pcounts <- lowpass <- numeric()
  pcpos <- lppos <- character()

  for (i in 1:NROW(output)){
    # relative abundance for each transcript
  	
    gpos <- rownames(etTable)%in%(v1names[p1==output[i,7]])
    et <- etTable[gpos,]
    #logcpm et[,2], min2[i] is the relative abundance. median + sd
    min2[i] <- median(et[is.finite(et[[2]]),2])+sd(et[is.finite(et[[2]]),2])
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
      if (fc[i]>2){ # if fold change > 2fold
#         if(length(maxpos[i]) != length(as.numeric(p2[v1names==rownames(passed[which.max(abs(passed[[1]])),])]))){
#           browser();
#         }
        maxpos[i] <- as.numeric(p2[v1names==rownames(passed[which.max(abs(passed[[1]])),])])[1]
        j <- maxpos[i]
        
        # get the relative position in transcript, 5UTR negative, CDS position, 3UTR positive + 100.
        if (j<=g[[2]]) {
          dist[i] <- j-g[[2]]-1
        } else if (j>=g[[3]]){
          dist[i] <- j-g[[3]]+1001
        } else dist[i] <- (j-g[[2]]+1)*1000/(g[[3]]-g[[2]]-1)
        passed_pos <- which(v1names%in%rownames(passed[abs(passed[[1]])>=1,])) ## only those positions that also pass FC filters
				
				#
        if (length(passed_pos)>1){
          passed_counts <- t(t(v1all[passed_pos,])*sf) # sf is from  yeast_combined_v1all_edgeR_sf.Rdata
          pcounts[i] <- as.numeric(length(which(rowSums(passed_counts>normalizedCountFilter)>=numConditionPassNCF))>0)
	        if (pcounts[i]>0) pcpos[i] <- paste(passed_pos[which(rowSums(passed_counts>normalizedCountFilter)>=numConditionPassNCF)],collapse=";")
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
	        if (pcounts[i]>0) pcpos[i] <- passed_pos;

          if (passed$logFC[abs(passed[[1]])>1]>0) { # increasing/decreasing of V1 signal 
            lowpass[i] <- as.numeric(min(passed_counts[group_names]) > mean(passed_counts[!group_names])+2*sd(passed_counts[!group_names]))
          } else {
            lowpass[i] <- as.numeric(max(passed_counts[group_names]) < mean(passed_counts[!group_names])-2*sd(passed_counts[!group_names]))            
          }
	        if (lowpass[i]>0) lppos[i] <- passed_pos
        }
      }
    }
  }

  output2 <- data.frame(output,maxPos=maxpos,maxFC=fc,m1sd=min2,dist=dist,pcounts,pcpos,lowpass,lppos)[fc>2,]
  print(paste(c,": ",NROW(output2[output2$pcounts>0&output2$lowpass>0,]),"positive",sep=""));
  if (NROW(output2[output2$pcounts>0&output2$lowpass>0,])>0){ 
     output2b <- output2[output2$pcounts>0&output2$lowpass>0,] 
     save(output2b,file=paste(resultdir,"combined_",c,"_output2_wfilters.Rdata",sep=""))
     write.table(output2b,file=paste(resultdir,"combined_",c,"_output2_wfilters.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
  }
}

  
