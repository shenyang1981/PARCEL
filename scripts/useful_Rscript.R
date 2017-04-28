readtab <- function(filename){
        v1 <- scan(filename,what="",sep="\n")
        v1 <- strsplit(v1, "\t")
        names(v1) <- sapply(v1, `[[`, 1)
        v1 <- lapply(v1, `[`, -1)
        v1 <- lapply(v1,function(x) as.numeric(strsplit(x,";")[[1]]))
        return(v1)}

scale_genic2 <- function(genic,l){
      n <- length(genic)
      if (1000%%n==0){
      	 return(rep(genic,rep(1000/n,n)))
      } else{
      	 g1k <- rep(genic,rep(l,n))/n
      	 steps <- seq(1,length(g1k),n)
      	 return(sapply(1:length(steps),function(i) sum(g1k[steps[i]:(steps[i]+n-1)])))}}

## sgl <- 1000
plotMetagene <- function(genes,probes,lmax){
             meta5 <- meta3 <- matrix(nrow=NROW(genes),ncol=500)
             metagene <- matrix(nrow=NROW(genes),ncol=sgl)
             for (i in 1:NROW(genes)){
                 if (genes[i,6]=="+"){
                    u5 <- (genes[i,2]-500):(genes[i,2]-1)
                    u3 <- (genes[i,3]+1):(genes[i,3]+500)
                    genic <- probes[genes[i,2]:genes[i,3]]
                    meta5[i,] <- c(rep(0,length(u5[u5<1])),probes[u5[u5>0&u5<=lmax]],rep(0,length(u5[u5>lmax])))
                    meta3[i,] <- c(rep(0,length(u3[u3<1])),probes[u3[u3>0&u3<=lmax]],rep(0,length(u3[u3>lmax])))
                 } else {
                    u5 <- (genes[i,3]+500+lmax):(genes[i,3]+1+lmax)
                    u3 <- (genes[i,2]-1+lmax):(genes[i,2]-500+lmax)
                    genic <- probes[(genes[i,3]+lmax):(genes[i,2]+lmax)]
                    meta5[i,] <- c(rep(0,length(u5[u5<=lmax])),probes[u5[u5>lmax&u5<=(lmax*2)]],rep(0,length(u5[u5>(lmax*2)])))
                    meta3[i,] <- c(rep(0,length(u3[u3<=lmax])),probes[u3[u3>lmax&u3<=(lmax*2)]],rep(0,length(u3[u3>(lmax*2)])))}
                 metagene[i,] <- scale_genic2(genic,sgl)}
             return(data.frame(meta5,metagene,meta3))}


plotMeta2 <- function(genes,probes,lmax){
          meta5 <- meta3 <- matrix(nrow=NROW(genes),ncol=1001)
          for (i in 1:NROW(genes)){
                 if (genes[i,6]=="+"){
                    u5 <- (genes[i,2]-500):(genes[i,2]+500)
                    u3 <- (genes[i,3]-500):(genes[i,3]+500)
                    meta5[i,] <- c(rep(0,length(u5[u5<1])),probes[u5[u5>0&u5<=lmax]],rep(0,length(u5[u5>lmax])))
                    meta3[i,] <- c(rep(0,length(u3[u3<1])),probes[u3[u3>0&u3<=lmax]],rep(0,length(u3[u3>lmax])))
                 } else {
                    u5 <- (genes[i,3]+500+lmax):(genes[i,3]-500+lmax)
                    u3 <- (genes[i,2]+500+lmax):(genes[i,2]-500+lmax)
                    meta5[i,] <- c(rep(0,length(u5[u5<=lmax])),probes[u5[u5>lmax&u5<=(lmax*2)]],rep(0,length(u5[u5>(lmax*2)])))
                    meta3[i,] <- c(rep(0,length(u3[u3<=lmax])),probes[u3[u3>lmax&u3<=(lmax*2)]],rep(0,length(u3[u3>(lmax*2)])))}}
          meta <- list()
          meta$meta5 <- meta5
          meta$meta3 <- meta3
          return(meta)}


kadane <- function(pos1,posN){
       max_temp <- max_now <- 0
       begin_temp <- begin <- end <- 0
       for (i in pos1:posN){
           if (max_temp < 0){
           max_temp <- score[i]
           begin_temp <- i} else {
              max_temp <- max_temp+score[i]}
           if (max_temp >= max_now){
              max_now <- max_temp
              begin <- begin_temp
              end <- i}}
        output <- cbind(begin=begin,end=end,maxScore=max_now)
        return (output)}

combine_list <- function(list){
             points <- NULL
             for (i in 1:NROW(list)){
                 points <- c(points,list[i,1]:list[i,2])}
                 return(points)}

annotateBac <- function(genes,output){
              gpos <- genes[genes[[6]]=="+",]
              gneg <- genes[genes[[6]]=="-",]

	      status <- geneID <- character()
	      dist <- numeric()
	      for (i in 1:NROW(output)){
    	      	  pos <- output[i,4]
    		  str <- output[i,5]
    		  status[i] <- "intergenic"
    		  geneID[i] <- dist[i] <- NA
    		  if (str=="+"){ # bedfile, start position - gpos[[2]]+1 
       		     d5 <- ifelse((gpos[[2]]+1)>pos,(gpos[[2]]+1)-pos,Inf)
       		     d3 <- ifelse(pos>gpos[[3]],pos-gpos[[3]],Inf)
		     if (min(d5)<500){
          	     	inter5 <- gpos[which.min(d5),]
          		geneID[i] <- as.character(inter5[,4])
          		status[i] <- "utr5"
          		dist[i] <- (inter5[,2]+1)-pos
       		     } else if (length(which(((gpos[[2]]+1)<=pos)&(gpos[[3]]>=pos)))>0){
          	        status[i] <- "genic"
         		ingene <- genes[which((genes[[6]]=="+")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)),]
          		geneID[i] <- as.character(ingene[1,4])
          		dist[i] <- pos-(ingene[1,2]+1)
       		     } else if (min(d3)<500){
          	        inter3 <- gpos[which.min(d3),]
          		geneID[i] <- as.character(inter3[,4])
          		status[i] <- "utr3"
         		dist[i] <- pos-inter3[,3]}
		  } else {
       		     d5 <- ifelse(gneg[[3]]<pos,pos-gneg[[3]],Inf)
       		     d3 <- ifelse(pos<(gneg[[2]]+1),(gneg[[2]]+1)-pos,Inf)
     		     if (min(d5)<500){
               	     	inter5 <- gneg[which.min(d5),]
          		geneID[i] <- as.character(inter5[,4])
          		status[i] <- "utr5"
          		dist[i] <- pos-inter5[,3]
       	 	     } else if (length(which(((gneg[[2]]+1)<=pos)&(gneg[[3]]>=pos)))>0){
          	        status[i] <- "genic"
          		ingene <- genes[which((genes[[6]]=="-")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)),]
          		geneID[i] <- as.character(ingene[1,4])
          		dist[i] <- ingene[1,3]-pos
       		     } else if (min(d3)<500){
          	        inter3 <- gneg[which.min(d3),]
          		geneID[i] <- as.character(inter3[,4])
          		status[i] <- "utr3"
          		dist[i] <- inter3[,2]+1-pos}}}
	      return(data.frame(output,status,geneID,dist))}

annotateBac_NC <- function(genes,output){
              gpos <- genes[genes[[6]]=="+",]
              gneg <- genes[genes[[6]]=="-",]
	      status2 <- geneID2 <- character()
	      dist2 <- numeric()
	      for (i in 1:NROW(output)){
    	      	  pos <- output[i,4]
    		  str <- output[i,5]
    		  status2[i] <- geneID2[i] <- dist2[i] <- NA
		  if (str=="+"){
       		     d5 <- ifelse((gpos[[2]]+1)>pos,(gpos[[2]]+1)-pos,Inf)
       		     d3 <- ifelse(pos>gpos[[3]],pos-gpos[[3]],Inf)
		     if (length(which(((gpos[[2]]+1)<=pos)&(gpos[[3]]>=pos)))>0){
		     	status2[i] <- "ncRNA"		     
			ingene <- genes[which((genes[[6]]=="+")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)),]	  
			geneID2[i] <- as.character(ingene[1,4])
			dist2[i] <- pos-(ingene[1,2]+1)}
    		  } else {
       		     d5 <- ifelse(gneg[[3]]<pos,pos-gneg[[3]],Inf)
       		     d3 <- ifelse(pos<(gneg[[2]]+1),(gneg[[2]]+1)-pos,Inf)  
       		     if (length(which(((gneg[[2]]+1)<=pos)&(gneg[[3]]>=pos)))>0){
          	     	status2[i] <- "ncRNA"
          		ingene <- genes[which((genes[[6]]=="-")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)),]
          		geneID2[i] <- as.character(ingene[1,4])
          		dist2[i] <- ingene[1,3]-pos}}}
    	      return(data.frame(output,status2,geneID2,dist2))}

assignGenes <- function(genes,etpassed){
	      points <- as.numeric(rownames(etpassed))
              gpos <- genes[genes[[6]]=="+",]
              gneg <- genes[genes[[6]]=="-",]

	      status <- geneID <- character()
	      dist <- numeric()
	      for (i in 1:length(points)){
	      	  str <- ifelse(points[i]>lmax,"-","+")
	      	  pos <- ifelse(str=="+",points[i],points[i]-lmax)
		  status[i] <- "intergenic"
    		  geneID[i] <- dist[i] <- NA
		  if (str=="+"){
		     d5 <- ifelse(gpos[[2]]>pos,gpos[[2]]-pos,Inf)
       		     d3 <- ifelse(pos>gpos[[3]],pos-gpos[[3]],Inf)
		     if (min(d5)<500){
          	     	inter5 <- gpos[which.min(d5),]
          		geneID[i] <- as.character(inter5[,4])
          		status[i] <- "utr5"
          		dist[i] <- inter5[,2]-pos
       	             } else if (length(which((genes[[6]]==str)&(genes[[2]]<=pos)&(genes[[3]]>=pos)))>0){
          	        status[i] <- "genic"
          		ingene <- genes[which((genes[[6]]==str)&(genes[[2]]<=pos)&(genes[[3]]>=pos)),]
			geneID[i] <- as.character(ingene[1,4])
          		dist[i] <- pos-ingene[1,2]
       		     } else if (min(d3)<500){
          	        inter3 <- gpos[which.min(d3),]
          		geneID[i] <- as.character(inter3[,4])
          		status[i] <- "utr3"
          		dist[i] <- pos-inter3[,3]}
    	     	  } else {
       		     d5 <- ifelse(gneg[[3]]<pos,pos-gneg[[3]],Inf)
       		     d3 <- ifelse(pos<gneg[[2]],gneg[[2]]-pos,Inf)
       		     if (min(d5)<500){
          	     	inter5 <- gneg[which.min(d5),]
          		geneID[i] <- as.character(inter5[,4])
          		status[i] <- "utr5"
          		dist[i] <- pos-inter5[,3]
       		     } else if (length(which((genes[[6]]==str)&(genes[[2]]<=pos)&(genes[[3]]>=pos)))>0){
          	        status[i] <- "genic"
          		ingene <- genes[which((genes[[6]]==str)&(genes[[2]]<=pos)&(genes[[3]]>=pos)),]
	  		geneID[i] <- as.character(ingene[1,4]) # if multiple, just take first one
          		dist[i] <- ingene[1,3]-pos
       		     } else if (min(d3)<500){
          	        inter3 <- gneg[which.min(d3),]
          		geneID[i] <- as.character(inter3[,4])
          		status[i] <- "utr3"
          		dist[i] <- inter3[,2]-pos}}}
return(data.frame(etpassed,points,status,geneID,dist))}

assignGenes2 <- function(genes,points){ ## directly from points
              # points <- as.numeric(rownames(etpassed))
              gpos <- genes[genes[[6]]=="+",]
              gneg <- genes[genes[[6]]=="-",]

              status <- geneID <- character()
              dist <- numeric()
              for (i in 1:length(points)){
                  str <- ifelse(points[i]>lmax,"-","+")
                  pos <- ifelse(str=="+",points[i],points[i]-lmax)
                  status[i] <- "intergenic"
                  geneID[i] <- dist[i] <- NA
                  if (str=="+"){
                     d5 <- ifelse(gpos[[2]]>pos,gpos[[2]]-pos,Inf)
                     d3 <- ifelse(pos>gpos[[3]],pos-gpos[[3]],Inf)
                     if (min(d5)<500){
                        inter5 <- gpos[which.min(d5),]
                        geneID[i] <- as.character(inter5[,4])
                        status[i] <- "utr5"
                        dist[i] <- inter5[,2]-pos
                     } else if (length(which((genes[[6]]==str)&(genes[[2]]<=pos)&(genes[[3]]>=pos)))>0){
                        status[i] <- "genic"
                        ingene <- genes[which((genes[[6]]==str)&(genes[[2]]<=pos)&(genes[[3]]>=pos)),]
			# if(NROW(ingene)>1) print(paste(pos,str,ingene[1,4],sep=":"))
                        geneID[i] <- as.character(ingene[1,4])
                        dist[i] <- pos-ingene[1,2]
                     } else if (min(d3)<500){
                        inter3 <- gpos[which.min(d3),]
                        geneID[i] <- as.character(inter3[,4])
                        status[i] <- "utr3"
                        dist[i] <- pos-inter3[,3]}
                  } else {
                     d5 <- ifelse(gneg[[3]]<pos,pos-gneg[[3]],Inf)
                     d3 <- ifelse(pos<gneg[[2]],gneg[[2]]-pos,Inf)
                     if(min(d5)<500){
                        inter5 <- gneg[which.min(d5),]
                        geneID[i] <- as.character(inter5[,4])
                        status[i] <- "utr5"
                        dist[i] <- pos-inter5[,3]
                     } else if(length(which((genes[[6]]==str)&(genes[[2]]<=pos)&(genes[[3]]>=pos)))>0){
                        status[i] <- "genic"
                        ingene <- genes[which((genes[[6]]==str)&(genes[[2]]<=pos)&(genes[[3]]>=pos)),]
			# if(NROW(ingene)>1) print(paste(pos,str,ingene[1,4],sep=":"))
                        geneID[i] <- as.character(ingene[1,4])
                        dist[i] <- ingene[1,3]-pos
                     } else if (min(d3)<500){
                        inter3 <- gneg[which.min(d3),]
                        geneID[i] <- as.character(inter3[,4])
                        status[i] <- "utr3"
                        dist[i] <- inter3[,2]-pos}}}
return(data.frame(points,geneID))}

## from positions, assign genes/status/dist, used in genomic_distribution scripts
## for utr5 and utr3, give positions, for genic, give rel. positions, gene length scaled to 1
assignG3 <- function(genes,points){ ## directly from points  
  gpos <- genes[genes[[6]]=="+",]
  gneg <- genes[genes[[6]]=="-",]
  status <- geneID <- character()
  dist <- numeric()
  for (i in 1:length(points)){
    str <- ifelse(points[i]>lmax,"-","+")
    pos <- ifelse(str=="+",points[i],points[i]-lmax)
    status[i] <- "intergenic"
    geneID[i] <- dist[i] <- NA
    if (str=="+"){
      d5 <- ifelse((gpos[[2]]+1)>pos,gpos[[2]]+1-pos,Inf)
      d3 <- ifelse(pos>gpos[[3]],pos-gpos[[3]],Inf)
      if (min(d5)<500){
        inter5 <- gpos[which.min(d5),]
        geneID[i] <- as.character(inter5[,4])
        status[i] <- "utr5"
        dist[i] <- pos-(inter5[,2]+1)
      } else if (length(which((genes[[6]]=="+")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)))>0){
        status[i] <- "genic"
        ingene <- genes[which((genes[[6]]=="+")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)),]
        if(NROW(ingene)>1) print(paste(pos,str,ingene[1,4],sep=":"))
        geneID[i] <- as.character(ingene[1,4])
        dist[i] <- (pos-(ingene[1,2]+1)+1)*1000/(ingene[1,3]-(ingene[1,2]+1)+1)
      } else if (min(d3)<500){
        inter3 <- gpos[which.min(d3),]
        geneID[i] <- as.character(inter3[,4])
        status[i] <- "utr3"
        dist[i] <- 1000+pos-inter3[,3]}      
    } else {
      d5 <- ifelse(gneg[[3]]<pos,pos-gneg[[3]],Inf)
      d3 <- ifelse(pos<(gneg[[2]]+1),gneg[[2]]+1-pos,Inf)
      if(min(d5)<500){
        inter5 <- gneg[which.min(d5),]
        geneID[i] <- as.character(inter5[,4])
        status[i] <- "utr5"
        dist[i] <- inter5[,3]-pos
      } else if(length(which((genes[[6]]=="-")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)))>0){
        status[i] <- "genic"
        ingene <- genes[which((genes[[6]]=="-")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)),]
        if(NROW(ingene)>1) print(paste(pos,str,ingene[1,4],sep=":"))
        geneID[i] <- as.character(ingene[1,4])
        dist[i] <- (ingene[1,3]-pos+1)*1000/(ingene[1,3]-(ingene[1,2]+1)+1)
      } else if (min(d3)<500){
        inter3 <- gneg[which.min(d3),]
        geneID[i] <- as.character(inter3[,4])
        status[i] <- "utr3"
        dist[i] <- 1000+(inter3[,2]+1)-pos}}}
  return(data.frame(points,geneID,status,dist))}

## assign positions to noncoding genes, if they are genic to NC genes
assignNC <- function(nc,points){   
  gpos <- genes[genes[[6]]=="+",]
  gneg <- genes[genes[[6]]=="-",]
  status <- geneID <- character()
  dist <- numeric()
  for (i in 1:length(points)){
    str <- ifelse(points[i]>lmax,"-","+")
    pos <- ifelse(str=="+",points[i],points[i]-lmax)
    status[i] <- "unassigned"
    geneID[i] <- dist[i] <- NA
    if (str=="+"){
      if (length(which((genes[[6]]=="+")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)))>0){
        status[i] <- "noncoding"
        ingene <- genes[which((genes[[6]]=="+")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)),]
        if(NROW(ingene)>1) print(paste(pos,str,ingene[1,4],sep=":"))
        geneID[i] <- as.character(ingene[1,4])
        dist[i] <- (pos-(ingene[1,2]+1)+1)*1000/(ingene[1,3]-(ingene[1,2]+1)+1)
      }} else {
    if(length(which((genes[[6]]=="-")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)))>0){
        status[i] <- "noncoding"
        ingene <- genes[which((genes[[6]]=="-")&((genes[[2]]+1)<=pos)&(genes[[3]]>=pos)),]
        if(NROW(ingene)>1) print(paste(pos,str,ingene[1,4],sep=":"))
        geneID[i] <- as.character(ingene[1,4])
        dist[i] <- (ingene[1,3]-pos+1)*1000/(ingene[1,3]-(ingene[1,2]+1)+1)}}}  
  return(data.frame(points,geneID,status,dist))}



morefilters <- function(output,genes,etTable,thr){
	    fc <- rpos <- numeric()
	    for (i in 1:NROW(output)){ ## get +/-500 from gene
    	    	w <- output[i,]
    		g <- genes[genes[[4]]==as.character(w[,9]),]
    		if (g[,6]=="+"){
       		   st <- g[,2]-500
		   ed <- g[,3]+500
       		} else {
       		   st <- g[,2]+lmax-500
		   ed <- g[,3]+lmax+500}

            	# positions in etTable in this gene
	      	et <- etTable[as.numeric(rownames(etTable))%in%(st:ed),]
	      	min2 <- median(et[is.finite(et[[2]]),2])+sd(et[is.finite(et[[2]]),2]) # relative abundance
	      	passed <- et[et[[2]]>=min2&et[[3]]<thr,]
	      	passed <- passed[as.numeric(rownames(passed))%in%c(w[,1]:w[,2]),]
	      	fc[i] <- 0	    
	      	if (NROW(passed)>0) fc[i] <-  2^(max(ifelse(passed[[1]]>0,passed[[1]],-passed[[1]])))
	      	if (fc[i]>2) rpos <- c(rpos,as.numeric(rownames(passed)))}   

	    out <- list()
	    out$output2 <- data.frame(output,maxFC=fc)[fc>2,]
	    out$rpos <- rpos
	    return(out)}

morefilters2 <- function(etPass,genes,etTable){
            min2 <- numeric()
	    for (i in 1:NROW(etPass)){ ## get +/-500 from gene
               g <- genes[genes[[4]]==as.character(etPass[i,6]),]
               if (g[,6]=="+"){
                   st <- g[,2]-500
                   ed <- g[,3]+500
                } else {
                   st <- g[,2]+lmax-500
                   ed <- g[,3]+lmax+500}
                # positions in etTable in this gene
                et <- etTable[as.numeric(rownames(etTable))%in%(st:ed),]
                min2[i] <- median(et[is.finite(et[[2]]),2])+sd(et[is.finite(et[[2]]),2])}	
	    filtered <- data.frame(etPass,min2)[abs(etPass[[1]])>1&etPass[[2]]>=min2,]
	    return(filtered)}			


checkStatus <- function(reg){
            ingene <- sapply(cdsReg,function(x) sum(as.numeric(reg%in%x)))
            if (length(which(ingene>0))>1) print(paste(i,"inside_multiple_genes",sep=":"))
            else if (length(which(ingene>0))==0) return("noncoding")
            else {
                 if(as.numeric(ingene[ingene>0])/length(reg)==1) return("coding")
                 else return("mixed")}}
genicExtend <- function(reg){
            if (length(reg)>=200) pos <- reg
            else {
                 ingene <- sapply(cdsReg,function(x) sum(as.numeric(reg%in%x)))
                 fgene <- cdsReg[which(ingene>0)][[1]] # found gene
                 mid <- floor(median(reg))
                 if (mid-99>=fgene[1]){
                    st <- mid-99
                    ed <- ifelse(mid+100<=fgene[length(fgene)],mid+100,fgene[length(fgene)])
                 } else {
                    st <- fgene[1]
                    ed <- ifelse(st+199<=fgene[length(fgene)],fgene[1]+199,fgene[length(fgene)])}
                 pos <- st:ed}
             if (pos[1]<=lmax){
                return(data.frame(start=pos[1],end=pos[length(pos)],strand="+"))
             } else return(data.frame(start=pos[1]-lmax,end=pos[length(pos)]-lmax,strand="-"))}


GXpt <- function(mid){
        ingene <- sapply(cdsReg,function(x) sum(as.numeric(mid%in%x)))
        fgene <- cdsReg[which(ingene>0)][[1]]
        if (mid-99>=fgene[1]){
                    st <- mid-99
                    ed <- ifelse(mid+100<=fgene[length(fgene)],mid+100,fgene[length(fgene)])
                 } else {
                    st <- fgene[1]
                    ed <- ifelse(st+199<=fgene[length(fgene)],fgene[1]+199,fgene[length(fgene)])}
        pos <- st:ed
        if (mid<=lmax){
                return(data.frame(start=pos[1],end=pos[length(pos)],strand="+"))
        } else return(data.frame(start=pos[1]-lmax,end=pos[length(pos)]-lmax,strand="-"))}


interExtend <- function(reg){
            lgene <- rgene <- NULL
            if (length(reg)>=200) pos <- reg
            else {
                 if (reg[1]<=lmax) {
                    left <- (reg[1]-cemp3)[reg[1]>cemp3]
                    if (length(left)>0) lgene <- cdsReg[reg[1]>cemp3][which.min(left)][[1]]
                    right <- (cemp2-reg[length(reg)])[cemp2<=lmax&cemp2>reg[length(reg)]]
                    if (length(right)>0) rgene <- cdsReg[cemp2<=lmax&cemp2>reg[length(reg)]][which.min(right)][[1]]
                 } else {
                    left <- (reg[1]-cemp3)[cemp3>lmax&reg[1]>cemp3]
                    if (length(left)>0) lgene <- cdsReg[cemp3>lmax&reg[1]>cemp3][which.min(left)][[1]]
                    right <- (cemp2-reg[length(reg)])[cemp2>reg[length(reg)]]
                    if (length(right)>0) rgene <- cdsReg[cemp2>reg[length(reg)]][which.min(right)][[1]]}
                 mid <- floor(median(reg))
                 if (is.null(lgene)) st <- ifelse(mid<=lmax,1,lmax+1)
                 else if (mid-99>lgene[length(lgene)]) st <- mid-99
                 else st <- lgene[length(lgene)]+1
                 if (is.null(rgene)) ed <- ifelse(mid<=lmax,lmax,lmax*2)
                 else ed <- ifelse(st+199<=rgene[1]-1,st+199,rgene[1]-1)
                 pos <- st:ed}
            if (pos[1]<=lmax){
                return(data.frame(start=pos[1],end=pos[length(pos)],strand="+"))
            } else return(data.frame(start=pos[1]-lmax,end=pos[length(pos)]-lmax,strand="-"))}

IXpt <- function(mid){
     lgene <- rgene <- NULL
     if (mid<=lmax) {
        left <- (mid-cemp3)[mid>cemp3]
        if (length(left)>0) lgene <- cdsReg[mid>cemp3][which.min(left)][[1]]
        right <- (cemp2-mid)[cemp2<=lmax&cemp2>mid]
        if (length(right)>0) rgene <- cdsReg[cemp2<=lmax&cemp2>mid][which.min(right)][[1]]
    } else {
        left <- (mid-cemp3)[cemp3>lmax&mid>cemp3]
        if (length(left)>0) lgene <- cdsReg[cemp3>lmax&mid>cemp3][which.min(left)][[1]]
        right <- (cemp2-mid)[cemp2>mid]
        if (length(right)>0) rgene <- cdsReg[cemp2>mid][which.min(right)][[1]]}
    if ((!is.null(lgene))&(!is.null(rgene))) {
       if (mid-99>lgene[length(lgene)]) st <- mid-99
       else st <- lgene[length(lgene)]+1
       ed <- ifelse(st+199<=rgene[1]-1,st+199,rgene[1]-1)
       pos <- st:ed
       if (pos[1]<=lmax){
       	  return(data.frame(start=pos[1],end=pos[length(pos)],strand="+"))
       } else return(data.frame(start=pos[1]-lmax,end=pos[length(pos)]-lmax,strand="-"))}}

splitCNC <- function(reg,flag){#flag=1,coding; flag=2,noncoding
         ingene <- sapply(cdsReg,function(x) sum(as.numeric(reg%in%x)))
         fgene <- cdsReg[which(ingene>0)][[1]]
         cpart <- reg[reg%in%fgene]
         if (flag==1) return(cpart)
         else {
              npart <- reg[!reg%in%fgene]
              n1 <- length(npart[npart<cpart[1]])
              n2 <- length(npart[npart>cpart[length(cpart)]])
              ## nparts spanning the cpart?
              if (n1>0&n2>0){
                 if (n1>n2)  return(npart[npart<cpart[1]])
                 else return(npart[npart>cpart[length(cpart)]])
              } else return(npart)}}


plot_coverage <- function(genes,pos){
              pts <- list()
              tss <- ifelse(genes[,6]=="+",genes[,2],genes[,3])
              tts <- ifelse(genes[,6]=="+",genes[,3],genes[,2])
              points <- pt5 <- pt3 <- numeric()
              for (i in 1:NROW(genes)){
                  if(length(pos[pos>=(genes[i,2]-500)&(pos<=genes[i,3]+500)])>0){
                        psd_pos <- pos[pos>=(genes[i,2]-500)&(pos<=genes[i,3]+500)]
                        if(length(psd_pos[psd_pos>=genes[i,2]&psd_pos<=genes[i,3]])>0){
                                psd_pos <- psd_pos[psd_pos>=genes[i,2]&psd_pos<=genes[i,3]]
                                points <- c(points,abs(psd_pos-tss[i]+1)*1000/(genes[i,3]-genes[i,2]+1))
                        } else if(genes[i,6]=="+"){
                          if(length(psd_pos[psd_pos<tss[i]])){
                                pt5 <- c(pt5,-tss[i]+psd_pos[psd_pos<tss[i]])
                          } else pt3 <- c(pt3,psd_pos[psd_pos>tts[i]]-tts[i]+1000)
                        } else {  # "-" stand
                          if(length(psd_pos[psd_pos>tss[i]])){
                                pt5 <- c(pt5,-psd_pos[psd_pos>tss[i]]+tss[i])
                          } else pt3 <- c(pt3,tts[i]-psd_pos[psd_pos<tts[i]]+1000)}}}
              pts$points <- points
              pts$pt5 <- pt5
              pts$pt3 <- pt3
              return(pts)}
