## adding filters to the called regions
# a). p-val < 10/n
# b). rel. abundance within transcript > median+sd
# c). abs(fold change) > 2
# d). normalized counts	 of two conditions > 10
# e). both test conditions pass mean+2sd of all control	conditions

setwd("/mnt/projects/sunm/others/yeast/analysis_combined")
source("/mnt/projects/sunm/bacteria/useful_Rscript.R")
load("/mnt/projects/sunm/bacteria/fastq/S288C4/sgdGenes_gStruc.Rdata")
v1 <- readtab("../fastq_combined/ala_1.tab")
l <- sapply(v1,length)+200
p1 <- as.character(unlist(mapply(rep,names(l),l)))
p2 <- as.character(unlist(lapply(l,function(x) c(1:x))))

load("yeast_combined_v1all.Rdata")
n <- length(which(rowSums(v1all)>NCOL(v1all)))
thr <- 10/n
v1names <- rownames(v1all)
load("yeast_combined_v1all_edgeR_sf.Rdata")

conditions <- c("arg","asn","asp","atp","cys","glucose","gly","his",
                "lys","oxal","pro","sam","ser","tpp","trp","tyr","utp","val")

for (c in conditions){
  print(c)
  load(paste("yeast_fastq2_",c,"_output10.Rdata",sep="")) # output
  load(paste("yeast_combined_v1all_etTable_",c,".Rdata",sep="")) # etTable
  group_names <- rep(FALSE,NCOL(v1all))
  group_names[grep(c,colnames(v1all))] <- TRUE
  
  fc <- maxpos <- dist <- min2 <- pcounts <- lowpass <- numeric()
  pcpos <- lppos <- character()

  for (i in 1:NROW(output)){
    gpos <- rownames(etTable)%in%(v1names[p1==output[i,7]])
    et <- etTable[gpos,]
    min2[i] <- median(et[is.finite(et[[2]]),2])+sd(et[is.finite(et[[2]]),2])
    passed <- et[et[[2]]>=min2[i]&et[[3]]<thr,]
    passed <- passed[rownames(passed)%in%v1names[output[i,1]:output[i,2]],]
    ## passed - pass min2, thr

    g <- gStruc[gStruc$gene==as.character(output[i,7]),]
    fc[i] <- maxpos[i] <- dist[i] <- pcounts[i] <- lowpass[i] <- 0  
    pcpos[i] <- lppos[i] <- "none"
    
    if (NROW(passed)>0){    
      fc[i] <- 2^(max(abs(passed[[1]])))
      if (fc[i]>2){
        maxpos[i] <- as.numeric(p2[v1names==rownames(passed[which.max(abs(passed[[1]])),])])
        j <- maxpos[i]
        if (j<=g[[2]]) {
          dist[i] <- j-g[[2]]-1
        } else if (j>=g[[3]]){
          dist[i] <- j-g[[3]]+1001
        } else dist[i] <- (j-g[[2]]+1)*1000/(g[[3]]-g[[2]]-1)
        passed_pos <- which(v1names%in%rownames(passed[abs(passed[[1]])>=1,])) ## only those positions that also pass FC filters

        if (length(passed_pos)>1){
          passed_counts <- t(t(v1all[passed_pos,])*sf)
          pcounts[i] <- as.numeric(length(which(rowSums(passed_counts>10)>1))>0)
	  if (pcounts[i]>0) pcpos[i] <- paste(passed_pos[which(rowSums(passed_counts>10)>1)],collapse=";")

	  passed2 <- passed[match(rownames(passed_counts),rownames(passed)),]  	     
	  sign <- ifelse(passed2[[1]]<0,TRUE,FALSE) ## TRUE for up, FALSE for down

	  lpstatus <- unlist(lapply(1:NROW(passed2),function(x){
		if (sign[x]){
		   return(min(passed_counts[x,group_names])>(mean(passed_counts[x,!group_names])+2*sd(passed_counts[x,!group_names])))
		} else return(max(passed_counts[x,group_names])<(mean(passed_counts[x,!group_names])-2*sd(passed_counts[x,!group_names])))
	  }))

	  lowpass[i] <- sum(as.numeric(lpstatus))>0
	  if (lowpass[i]>0) lppos[i] <- paste(passed_pos[lpstatus],collapse=";")

        } else {
          passed_counts <- v1all[passed_pos,]*sf
          pcounts[i] <- as.numeric(sum(passed_counts>10)>1)
	  if (pcounts[i]>0) pcpos[i] <- passed_pos

          if (passed$logFC[abs(passed[[1]])>=1]<0) { #


# sign <- ifelse(passed2[[1]]>0,TRUE,FALSE)
(passed$logFC[abs(passed[[1]])>1]>0)





            lowpass[i] <- as.numeric(min(passed_counts[group_names]) > mean(passed_counts[!group_names])+2*sd(passed_counts[!group_names]))
          } else {
            lowpass[i] <- as.numeric(max(passed_counts[group_names]) < mean(passed_counts[!group_names])-2*sd(passed_counts[!group_names]))            
          }
	  if (lowpass[i]>0) lppos[i] <- passed_pos}
}}}

    output2 <- data.frame(output,maxPos=maxpos,maxFC=fc,m1sd=min2,dist=dist,pcounts,pcpos,lowpass,lppos)[fc>2,]
    print(NROW(output2[output2$pcounts>0&output2$lowpass>0,]))
    if (NROW(output2[output2$pcounts>0&output2$lowpass>0,])>0){ 
       output2b <- output2[output2$pcounts>0&output2$lowpass>0,] 
       save(output2b,file=paste("yeast_combined_",c,"_output2_wfilters.Rdata",sep=""))
}}

  
