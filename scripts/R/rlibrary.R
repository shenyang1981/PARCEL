suppressWarnings(suppressPackageStartupMessages(library(adagio)));

getWeightFold = function(siteinfo=NULL,siginfo=NULL,genelocus=NULL,geneinfo=NULL,positiveRegion=NULL,ext.win.size=100){
  genelocus.sort = bedr.sort.region(genelocus,check.chr=F,verbose=F,check.merge=F);
  siteinfo.sort = bedr.sort.region(siteinfo,check.chr=F,verbose=F,check.merge=F);
  siginfo.sort = bedr.sort.region(siginfo,check.chr=F,verbose=F,check.merge=F);
  
  # get region closest gene
  regioninGene = bedr(
    engine = "bedtools", 
    input = list(a = siginfo.sort,b=genelocus.sort),
    method = "closest", 
    params = "-s -D b",
    check.chr=F,
    check.merge=F,
    verbose=F
  );
  #browser();
  regioninGene = regioninGene[,c(1,2,3,4,5,6,10,13)]
  
  colnames(regioninGene) = c("chr","start","end","region","evalue","strand","gene","distance");
  
  # drop region without flanked gene
  regioninGene = regioninGene[regioninGene[,"gene"]!=".",];
  
  regioninGene[,"evalue"] = as.numeric(regioninGene[,"evalue"]);
  
  # all.x = T means keep region without flanked genes
  #regioninGeneFold = merge(regioninGene,geneinfo[,c("chr","gene","logfold")],by=c("gene","chr"),all.x=T);
  
  # all.x = F will drop region without flanked genes
  regioninGeneFold = merge(regioninGene,geneinfo[,c("chr","gene","logfold")],by=c("gene","chr"),all.x=F);
  
  
  siginfo.sort.extend = siginfo.sort;
  sigMidpoint = round((siginfo.sort.extend[,2] + siginfo.sort.extend[,3])/2,0);
  sigWidth = siginfo.sort.extend[,3] - siginfo.sort.extend[,2] + 1;
  siginfo.sort.extend[sigWidth<ext.win.size*2,2] = ifelse(sigMidpoint[sigWidth<ext.win.size*2] -ext.win.size>0,sigMidpoint[sigWidth<ext.win.size*2] -ext.win.size,1);
  siginfo.sort.extend[sigWidth<ext.win.size*2,3] = sigMidpoint[sigWidth<ext.win.size*2] +ext.win.size;
  colnames(siginfo.sort.extend)[1:5] = c("chr","start","end","V4","V5");
  
  siginfo.sort.extend = bedr.sort.region(siginfo.sort.extend,check.chr=F,verbose=F,check.merge=F);
  
  siteinSigRegion = bedr(
    engine = "bedtools", 
    input = list(a = siginfo.sort.extend ,b=siteinfo.sort), 
    method = "intersect", 
    params = "-s -loj -sorted",
    check.chr=F,
    verbose=F
  );
# 
# siteinSigRegion = siteinSigRegion[,c("V1","V2","V3","V4","V5","V9","V10","V11")];
  
  #siteinSigRegion = bedr.join.region(siginfo.sort.extend,siteinfo.sort,check.chr=F,verbose=F);
  #browser();
  siteinSigRegion = siteinSigRegion[,c("chr","start","end","V4","V5","strand","P","logfold","sumcount")];

  colnames(siteinSigRegion) = c("chr","start","end","region","evalue","strand","pvalue","sitefold","sumcount");
  siteinSigRegion[,"evalue"] = as.numeric(siteinSigRegion[,"evalue"]);
  siteinSigRegion[,"pvalue"] = as.numeric(siteinSigRegion[,"pvalue"]);
  siteinSigRegion[,"sumcount"] = as.numeric(siteinSigRegion[,"sumcount"]);
  siteinSigRegion[,"sitefold"] = as.numeric(siteinSigRegion[,"sitefold"]);
  
  #browser();
  #allInfo = merge(siteinSigRegion[,c("region","strand","pvalue","sitefold","sumcount")],regioninGeneFold[,c("region","strand","gene","logfold")],by=c("region","strand"),all.x=T);
  
  allInfo = merge(siteinSigRegion[,c("region","strand","pvalue","sitefold","sumcount")],regioninGeneFold[,c("region","strand","gene","logfold")],by=c("region","strand"),all.x=F);

  weightFold = function(x){
    pcut = 1e-1;
    #weight = x[,"sumcount"]/sum(x[,"sumcount"]);
    weight = log(x[,"sumcount"]+1)/sum(log(x[,"sumcount"]+1));
    correctFold = x[,"sitefold"]-x[,"logfold"]; # correct site fold by subtracting site fold from gene fold
    if(sum(x[,"pvalue"]<pcut)>1){
      y = sum(abs(correctFold[x[,"pvalue"]<pcut]) * weight[x[,"pvalue"]<pcut]);
      ymax = max(abs(correctFold[x[,"pvalue"]<pcut]));
      #y = mean(abs(correctFold[x[,"pvalue"]<0.05]));
      #y = max(abs(correctFold[x[,"pvalue"]<0.05]));
    } else {
      y=0;
      ymax=0;
    }
    result = c(y,ymax)
    return(result);
  }
  myweightFold = ddply(allInfo, c("region"), weightFold);
  colnames(myweightFold) = c("region","weightFold","maxFold");
  finallist = merge(myweightFold,regioninGeneFold,by="region",all.x=T);
  if(!is.null(positiveRegion)){
    positive.sort = bedr.sort.region(positiveRegion,check.chr=F,verbose=F,check.merge=F);
    regioninPos = bedr.join.region(siginfo.sort.extend,positive.sort,check.chr=F,verbose=F,check.merge=F);
    #browser();
    regioninPos = regioninPos[,c(4,12)]
    colnames(regioninPos) = c("region","posRegion");  
    finallist = merge(finallist,regioninPos,by="region",all.x=T);
  }
  #browser();
  finallist[,"logE"] = log10(finallist[,"evalue"])*-1;
  finallist[finallist[,"logE"]>30,"logE"] = 30;
  return(finallist);
}


compareBatch = function(mycountData=NULL,coninfo=NULL,cond1=NULL,cond2=NULL,diffmethod="DEseq",prefix="",mincov=20,bmincov=40){
  
  realdeseq = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = mincov,method = diffmethod,coninfo=coninfo,autobatch = T,bmincov=bmincov);
  
  realdeseqNoB = diffSite(mycount = mycountData,cond1 = cond1,cond2 = cond2,maxrsum = mincov,method = diffmethod,coninfo=coninfo,autobatch = F,bmincov=bmincov);
  
  all = merge(realdeseq,realdeseqNoB,by="id");
  all[,"id"] = as.numeric(as.character(all[,"id"]));
  all[,"start"] = all[,"id"];
  all[,"end"] = all[,"id"]+1;
  all[,"chr"] = "chr";
  all[,"locus"] = bedToLocus(all[,c("chr","start","end")]);
  
  posR = read.table("sumResult/positiveSAM.bed",header=F,stringsAsFactors = F);
  posR[,2] = posR[,2]-100;
  posR[,3] = posR[,3]+100;
  
  png(paste(prefix,"_Manhattan.png",sep=""),width = 1200,height = 800,pointsize=20,res=72);
  plot(all[,"id"],all[,"pvalue.x"],
       ylim=range(c(all[,"pvalue.x"],all[,"pvalue.y"])),
       cex=0.2,col="grey",xlab="genomic coordinate",ylab="-log10(P)");
  abline(h=3,lty=2,lwd=0.3,col="grey");
  abline(v=posR[!grepl("Not",posR[,4]),2],col="grey",lwd=0.5)
  abline(v=posR[grepl("Not",posR[,4]),2],col=trop[4],lwd=0.5,lty=2)
  points(all[all[,"pvalue.x"]>3,"id"],all[all[,"pvalue.x"]>3,"pvalue.x"],cex=0.3,col=trop[1])
  points(all[all[,"pvalue.y"]>3,"id"],all[all[,"pvalue.y"]>3,"pvalue.y"],cex=0.3,col=trop[2])
  legend("topright",legend=c("Considering Batch","Not Considering Batch"),pch=c(20,20),col=c(trop[1],trop[2]));
  dev.off();
  
  png(paste(prefix,"_Manhattan_Batch.png",sep=""),width = 1200,height = 800,pointsize=20,res=72);
  plot(all[,"id"],all[,"pvalue.x"],
       ylim=range(c(all[,"pvalue.x"],all[,"pvalue.y"])),
       cex=0.2,col="grey",xlab="genomic coordinate",ylab="-log10(P)");
  abline(h=3,lty=2,lwd=0.3,col="grey");
  abline(v=posR[!grepl("Not",posR[,4]),2],col="grey",lwd=0.5)
  abline(v=posR[grepl("Not",posR[,4]),2],col=trop[4],lwd=0.5,lty=2)
  points(all[all[,"pvalue.x"]>3,"id"],all[all[,"pvalue.x"]>3,"pvalue.x"],cex=0.3,col=trop[1])
  dev.off();
  
  png(paste(prefix,"_Manhattan_NotBatch.png",sep=""),width = 1200,height = 800,pointsize=20,res=72);
  plot(all[,"id"],all[,"pvalue.x"],
       ylim=range(c(all[,"pvalue.x"],all[,"pvalue.y"])),
       cex=0.2,col="grey",xlab="genomic coordinate",ylab="-log10(P)");
  abline(h=3,lty=2,lwd=0.3,col="grey");
  abline(v=posR[!grepl("Not",posR[,4]),2],col="grey",lwd=0.5)
  abline(v=posR[grepl("Not",posR[,4]),2],col=trop[4],lwd=0.5,lty=2)
  points(all[all[,"pvalue.y"]>3,"id"],all[all[,"pvalue.y"]>3,"pvalue.y"],cex=0.3,col=trop[2])
  dev.off();
  
  posR.sort = bedr.sort.region(posR[!grepl("Not",posR[,4]),c(1:4)],check.chr=F);
  result.sort = bedr.sort.region(all[,c("chr","start","end")],check.chr=F);
  is.region <- in.region(result.sort, posR.sort,check.chr=F);
  overlapSite = result.sort[is.region,]
  
  allPos = all[!is.na(match(all[,"locus"],bedToLocus(overlapSite))),]
  plot(allPos[,3],allPos[,6],xlab="-log10(P) Considering Batch",ylab="-log10(P) Not Considering Batch",pch=20,col=trop[1],xlim=range(c(allPos[,3],allPos[,6])),ylim=range(c(allPos[,3],allPos[,6])),main="Sites within Known RiboSwitch Region",cex=0.7);
  abline(0,1,col=trop[4],lty=2);
  
  allNonePos = all[is.na(match(all[,"locus"],bedToLocus(overlapSite))),]
  smoothScatter(allNonePos[,3],allNonePos[,6],xlab="-log10(P) Considering Batch",ylab="-log10(P) Not Considering Batch",pch=20,col=trop[2],xlim=range(c(allNonePos[,3],allNonePos[,6])),ylim=range(c(allNonePos[,3],allNonePos[,6])),main="Sites Not in Known RiboSwitch Region",cex=0.5);
  abline(0,1,col=trop[4],lty=2);
  return(all);
}


solve.default = function (a, b, tol = 1e-30, LINPACK = FALSE, ...) 
{
  if (is.complex(a) || (!missing(b) && is.complex(b))) {
    a <- as.matrix(a)
    if (missing(b)) {
      b <- diag(1 + (0+0i), nrow(a))
      colnames(b) <- rownames(a)
    }
    return(.Internal(La_solve_cmplx(a, b)))
  }
  if (is.qr(a)) {
    warning("solve.default called with a \"qr\" object: use 'qr.solve'")
    return(solve.qr(a, b, tol))
  }
  a <- as.matrix(a)
  if (missing(b)) {
    b <- diag(1, nrow(a))
    colnames(b) <- rownames(a)
  }
  .Internal(La_solve(a, b, tol))
}

unlockBinding("solve.default", as.environment("package:base"))
assign("solve.default", solve.default, "package:base")

v1s1ratio = function(mycount=NULL,cond1=NULL,cond2=NULL,coninfo=NULL,isexactTest=F,batcheffect=F,maxrsum=4,method="DEseq",nthreads=1,autobatch=F,bmincov=1){
  
  suppressPackageStartupMessages(library(edgeR));
  suppressPackageStartupMessages(library(DESeq2));
  suppressPackageStartupMessages(library(sva));
  # mycount is count matrix
  # cond1 vs. cond2, cond1 is control, cond2 is treatment
  # coninfo is data.frame with two columns: "cond" and "cols"
  # select rows with sum > maxrsum
  coninfo = coninfo[!is.na(match(coninfo[,"cond"],c(cond1,cond2))),];
  coninfo[,"cols"] = as.character(coninfo[,"cols"]);
  rownames(coninfo) = coninfo[,"cols"];
  mycountData = mycount[,coninfo[,"cols"]];
  mycountData = mycountData[rowSums(mycountData)>maxrsum,];
  #   colData  = data.frame("cond" = coninfo[,"cond"]);
  #   rownames(colData) = coninfo[,"cols"];
  colData = coninfo;
  finalresult = NULL;
  #if(autobatch==T){
    cond = coninfo[,"cond"];
    #     dds = DESeqDataSetFromMatrix(
    #       countData = mycountData[rowVars(mycountData)>0,],
    #       colData=colData,
    #       design = ~ cond);
    #     dds <- estimateSizeFactors(dds);
    #     dat <- counts(dds, normalized=TRUE)
    #     mod1 <- model.matrix(~ cond, colData(dds))
    #     mod0 <- model.matrix(~ 1, colData(dds))
    #     batch_unsup_sva = svaseq(dat, mod1, mod0);
    dds <- DGEList(counts=mycountData, group=cond);
    dds <- calcNormFactors(dds, method="upperquartile")
    dat = cpm(dds,log =F,prior.count=0);
    rawcount = dds$counts;
    rawcount = rawcount[rowSums(dat)>=1,];
    dat = dat[rowSums(dat)>=1,]; # 1 per million

    #dat = apply(dat,2,function(x){y = log2(x+1)});
    #finalresult = data.frame("id" = rownames(etTable), "sumcount" = rowSums(mycountData),"pvalue"= -1*log10(etTable[,"PValue"]),"logfold"=etTable[,"logFC"]);  
  #}
  #browser();
  cond1samples = coninfo[!is.na(match(coninfo[,"cond"],cond1)),"cols"];
  cond2samples = coninfo[!is.na(match(coninfo[,"cond"],cond2)),"cols"];
  v1s1 = apply(dat[,c(cond1samples,cond2samples)],1,function(x){
    y1 = log2(mean(x[cond1samples])+1);
    y2 = log2(mean(x[cond2samples])+1);
    result = y1 - y2
    return(result);
  });
  finalresult = data.frame("id" = rownames(dat),"ratio"=v1s1);
  myfinalresult = list("result"=finalresult,"norcpm"=dat,"rawcount"=rawcount);
  return(myfinalresult);
}


diffSite = function(mycount=NULL,cond1=NULL,cond2=NULL,coninfo=NULL,isexactTest=F,batcheffect=F,maxrsum=20,method="DEseq",nthreads=1,autobatch=F,bmincov=1){
  
  suppressPackageStartupMessages(library(edgeR));
  suppressPackageStartupMessages(library(DESeq2));
  suppressPackageStartupMessages(library(sva));
  # mycount is count matrix
  # cond1 vs. cond2, cond1 is control, cond2 is treatment
  # coninfo is data.frame with two columns: "cond" and "cols"
  # select rows with sum > maxrsum
  coninfo = coninfo[!is.na(match(coninfo[,"cond"],c(cond1,cond2))),];
  coninfo[,"cols"] = as.character(coninfo[,"cols"]);
  rownames(coninfo) = coninfo[,"cols"];
  mycountData = mycount[,coninfo[,"cols"]];
  mycountData = mycountData[rowSums(mycountData)>maxrsum,];
  #   colData  = data.frame("cond" = coninfo[,"cond"]);
  #   rownames(colData) = coninfo[,"cols"];
  colData = coninfo;
  finalresult = NULL;
  if(autobatch==T){
    cond = coninfo[,"cond"];
#     dds = DESeqDataSetFromMatrix(
#       countData = mycountData[rowVars(mycountData)>0,],
#       colData=colData,
#       design = ~ cond);
#     dds <- estimateSizeFactors(dds);
#     dat <- counts(dds, normalized=TRUE)
#     mod1 <- model.matrix(~ cond, colData(dds))
#     mod0 <- model.matrix(~ 1, colData(dds))
#     batch_unsup_sva = svaseq(dat, mod1, mod0);
    dds <- DGEList(counts=mycountData, group=cond);
    #dat <- calcNormFactors(dds, method="upperquartile")$counts
    dds <- calcNormFactors(dds, method="upperquartile")
    dat = cpm(dds);
    mod1 <- model.matrix(~ cond)
    mod0 <- cbind(mod1[,1])
    #browser();
    batch_unsup_sva = svaseq(dat[rowVars(dat)>0 & rowSums(dat)>bmincov,], mod1, mod0);
    if(batch_unsup_sva$n.sv!=0){
      # get normalized data with batch effect removed
#       design <- model.matrix(~1,data=cond)
#       combatNor <- ComBat(dat= logCPM,batch=batch_unsup_sva$sv,mod=design,par.prior=TRUE);

      if(batch_unsup_sva$n.sv>1){ # number of surrogate variable > 1
        for(i in 1:batch_unsup_sva$n.sv){
          coninfo[,paste("batch",i,sep="")] = batch_unsup_sva$sv[,i];  
          stripchart(batch_unsup_sva$sv[,i] ~ colnames(mycountData),vertical=TRUE,main=paste("SV",i,sep=""));
          abline(h=0)
        }  
      } else {
        #plot(batch_unsup_sva$sv,pch=19,main="unsupervised sva")
        stripchart(batch_unsup_sva$sv ~ colnames(mycountData),vertical=TRUE,main="SV1")
        abline(h=0)
        coninfo[,"batch"] = batch_unsup_sva$sv;  
      }
      colData = coninfo;
      batcheffect=T
    } else{
      batcheffect=F;
    }
  }
  if(method == "DEseq"){
    if(nthreads>1){
      suppressPackageStartupMessages(library("BiocParallel"));
      register(MulticoreParam(nthreads))
    }
    if(batcheffect==T){
      dds <- DESeqDataSetFromMatrix(countData = mycountData,
                                    colData = colData,
                                    design = ~ cond
                                    )
      tmp = colnames(colData);
      batchcols = tmp[grepl("batch",tmp)];
      design(dds) = as.formula(paste("~","cond + ",paste(batchcols,collapse=" + "),sep= " "));
      #browser();
    } else {
      dds <- DESeqDataSetFromMatrix(countData = mycountData,
                                    colData = colData,
                                    design = ~ cond)
    }
    logCPM=assay(rlog(dds));
    dds <- DESeq(dds)
    result <- results(dds, contrast=c("cond",cond1,cond2),cooksCutoff=F) ## contrast specifies conditions to be tested
    finalresult = data.frame("id" = rownames(result),"sumcount" = rowSums(mycountData),"pvalue" = -1*log10(result[,"pvalue"]),"logfold" = result[,"log2FoldChange"],"baseExp"= result[,"baseMean"]);
  }
  if(method == "edgeR"){
    if(batcheffect==T){
      # cond = coninfo[,"cond"];
      # batch = coninfo[,"batch"];
      #design <- model.matrix(~ batch + cond);
            
      tmp = colnames(colData);
      batchcols = tmp[grepl("batch",tmp)];
      myformula= as.formula(paste("~",paste(batchcols,collapse=" + ")," + cond",sep= " "));
      design <- model.matrix(data = colData,myformula);
      
      
      rownames(design) <- coninfo[,"cols"];
      y <- DGEList(counts = mycountData, genes = rownames(mycountData));
      y <- calcNormFactors(y)
      logCPM <- cpm(y, log=TRUE, prior.count=0.1); # normalized data
      y = estimateDisp(y, design,robust=T)
      fit <- glmQLFit(y, design, robust=TRUE)
      lrt = glmQLFTest(fit);
      etTable = lrt$table;
    }else{
      if(isexactTest==T){
        y = DGEList(counts = mycountData,group=coninfo[,"cond"],genes=rownames(mycountData));
        y <- calcNormFactors(y);
        logCPM <- cpm(y, log=TRUE, prior.count=0.1); # normalized data
        #y = estimateDisp(y,robust=T);
        y <- estimateCommonDisp(y, verbose = TRUE)
        y <- estimateTagwiseDisp(y)
        
        et <- exactTest(y, pair = c(cond2,cond1), dispersion = "auto")
        etTable <- et$table
      } else {
        cond = coninfo[,"cond"];
        design <- model.matrix(~ cond);
        rownames(design) <- coninfo[,"cols"];
        #y = DGEList(counts = mycountData,group=coninfo[,"cond"],genes=rownames(mycountData));
        y = DGEList(counts = mycountData,genes=rownames(mycountData));
        y <- calcNormFactors(y)
        logCPM <- cpm(y, log=TRUE, prior.count=0.1); # normalized data
        #y = estimateDisp(y,design,robust=T);
        
        y <- estimateGLMCommonDisp(y,design)
        y <- estimateGLMTrendedDisp(y,design)
        y <- estimateGLMTagwiseDisp(y,design)

        fit <- glmQLFit(y, design, robust=TRUE);
        lrt = glmQLFTest(fit);
        etTable = lrt$table;
      }
    }
    finalresult = data.frame("id" = rownames(etTable), "sumcount" = rowSums(mycountData),"pvalue"= -1*log10(etTable[,"PValue"]),"logfold"=etTable[,"logFC"],"baseExp"=etTable[,"logCPM"]);  
  }
  #browser();
  cond1samples = coninfo[!is.na(match(coninfo[,"cond"],cond1)),"cols"];
  cond2samples = coninfo[!is.na(match(coninfo[,"cond"],cond2)),"cols"];
  minchange = apply(logCPM[,c(cond1samples,cond2samples)],1,function(x){
    y1 = mean(x[cond1samples]);
    y2 = x[cond2samples] - y1;
    samedirection = abs(sum(y2)) == sum(abs(y2));
    result = ifelse(samedirection==T,min(abs(y2)),0);
    return(result);
    });
  
  myfinalresult = list("result"=finalresult,"norcpm"=logCPM,"rawcount"=mycountData,"minchange"=minchange);
  return(myfinalresult);
}

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

kadaneSun = function(score=NULL){
  pos1 <- 1
  posN <- length(score)
  output <- kadane(pos1=pos1,posN=posN)
  outtemp <- cbind(output,label=1)
  count <- NROW(output)
  while (output[NROW(output),3] > 5){
    left <- kadane(pos1=pos1,posN=(output[NROW(output),1]-1))
    right <- kadane(pos1=(output[NROW(output),2]+1),posN=posN)
    outtemp <- rbind(outtemp,cbind(left,label=0),cbind(right,label=0))
    outtemp <- outtemp[order(outtemp[,1]),]  
    maxBegin <- outtemp[outtemp[,4]==0,1][which.max(outtemp[outtemp[,4]==0,3])]
    maxInd <- which(outtemp[,1]==maxBegin)
    output <- rbind(output,outtemp[maxInd,1:3])
    pos1 <- ifelse(maxInd>1,outtemp[maxInd-1,2]+1,1)
    if (maxInd<NROW(outtemp)){
      posN <- outtemp[maxInd+1,1]-1
    } else posN <- length(score)
    outtemp[maxInd,4] <- 1
    count <- count+1
  }
  return(output);
}

kadaneShenSlow = function(x,pos=NULL,minscore=0,dscore=0){
  if(!is.null(pos) & length(x) == length(pos)){
    tmp = rep(dscore,max(pos));
    tmp[pos]=x;
  } else {
    tmp = x;
  }
  noarray = F;
  result = NULL;
  while(!noarray){
    tmpsub = maxsub(tmp,inds=T,compiled=T);
    if(tmpsub$sum < minscore){
      noarray=T;
    } else {
      if(is.null(result)){
        result = data.frame("begin"=tmpsub$inds[1],"end"=tmpsub$inds[2],"maxScore"=tmpsub$sum);
      } else {
        result = rbind(result,data.frame("begin"=tmpsub$inds[1],"end"=tmpsub$inds[2],"maxScore"=tmpsub$sum));
      }
      tmp[tmpsub$inds[1]:tmpsub$inds[2]]=dscore;
    }
  }
  return(result);
}

kadaneShen = function(x,pos=NULL,minscore=0,dscore=0,poslist=NULL){
  if(!is.null(pos) & length(x) == length(pos)){
    tmp = rep(dscore,max(pos)-min(pos)+1);
    tmp[pos-min(pos)+1]=x;
    pos = min(pos):max(pos);
  } else {
    tmp = x;
    pos=1:length(x);
  }
  tmpsub = maxsub(tmp,inds=T,compiled=T);
  if(tmpsub$sum < minscore){
      if(is.null(poslist)){
        poslist[[1]] = list("start"=1,"end"=2,"maxScore"=-10);
      }
      return(poslist);
  } else {
     # if there are maximum subarray, split array into two sub array
    leftpos = tmpsub$inds[1]+min(pos)-1;
    rightpos = tmpsub$inds[2]+min(pos)-1;
    resultlist = list("start"=leftpos,"end"=rightpos,"maxScore"=tmpsub$sum);
   
    # store current max window
    if(is.null(poslist)){
        poslist=list();
        poslist[[1]] = resultlist;
      } else {
       poslist[[length(poslist)+1]]=resultlist;
    }
    
     # handle left subarray
    if(tmpsub$inds[1]>1){
        xleft = tmp[1:(tmpsub$inds[1]-1)];
        posleft = pos[1:(tmpsub$inds[1]-1)]; # real position of left subarray
        if(length(xleft)==1 & sum(xleft)>minscore){
          poslist[[length(poslist)+1]] = list("start"=posleft[1],"end"=posleft[1],"maxScore"=xleft);
        } else {
          poslist = kadaneShen(xleft,pos=posleft,minscore=minscore,dscore=dscore,poslist=poslist);    
        }
      }
          
     # handle right subarray
    if(tmpsub$inds[2]<length(tmp)){
          xright = tmp[(tmpsub$inds[2]+1):length(tmp)];
          posright = pos[(tmpsub$inds[2]+1):length(tmp)]; # real position of left subarray
          if(length(xright)==1 & sum(xright)>minscore){
            poslist[[length(poslist)+1]] = list("start"=posright[1],"end"=posright[1],"maxScore"=xright);
          } else {
            poslist = kadaneShen(xright,pos=posright,minscore=minscore,dscore=dscore,poslist=poslist);    
          }
        }
      
  }
  return(poslist);
}

evalueOld = function(maxScoreBin=NULL,diffSiteResult = NULL,mincov=6){
  lambda = 0.862871
  K = 0.0809635
  n = sum(diffSiteResult[,"sumcount"]>mincov);
  tmp = diffSiteResult[diffSiteResult[,"sumcount"]>mincov,];
  maxScoreBin[,"Evalue"] = apply(maxScoreBin,1,function(x){
    y1 = tmp[tmp[,"pos"]>=x[1] & tmp[,"pos"]<=x[2],"score"];
    X <- sum(y1)-log(n,base=exp(1))/lambda;
    Prob <- K*exp(-lambda*X);
    return(Prob);
  })
  return(maxScoreBin);
}

evalue = function(n = NULL,sumscore=NULL){
  lambda = 0.862871
  K = 0.0809635
  X <- sumscore-log(n,base=exp(1))/lambda;
  Prob <- K*exp(-lambda*X);
  return(Prob);
}

test = function(x,mylist=NULL){
  if(x>10){
    return(mylist);
  } else {
    if(is.null(mylist)){
      mylist[0]=x;
    } else {
      mylist[length(mylist)+1]=x;  
    }
    test(x+1,mylist=mylist);
  }
}

percentGenomeCov = function(x = NULL,genomesize = NULL, covcutoff = NULL, bin=NULL , coordinates = NULL,maxcov=F){
  # given coverage cut off, how many genome is covered. 
  #x = covinfoW[,3];
  if(!is.null(bin) & !is.null(coordinates)){
    mybin = floor(coordinates/bin);
    data_t = data.table("bin" = mybin, "cov"=x);
    if(maxcov==T){
      countPerBin = data_t[,list("total"=max(cov)),by="bin"];  
    } else {
      countPerBin = data_t[,list("total"=sum(cov)),by="bin"];
    }
    tmp = as.data.frame(table(countPerBin$total));
    genomesize = floor(genomesize/bin);
  } else {
    tmp = as.data.frame(table(x))  
  }
  tmp[,1] = as.numeric(as.character(tmp[,1])); 
  tmp = tmp[order(tmp[,1]),];
  tmp[1,2] = tmp[1,2] + genomesize - sum(tmp[,2]); # zero coverage
  tmp[,3] = tmp[,2]/genomesize;
  tmp[,4] = cumsum(tmp[,3]);
  y = 1 - tmp[tmp[,1]>=covcutoff,][1,4];
  return(y);
}

slidingCov = function(x=NULL,position=NULL,windowsize = 300, stepsize = 50){
  tmp = array()
  tmp[position]=x;
  tmp[is.na(tmp)] = 0;
  result = rollapply(tmp,width = windowsize,by = stepsize,mean);
  return(result);
}

downSample = function(x=NULL,prop = NULL){
  y = rbinom(n = length(x),size=x,p=prop);
  return(y);
}


simulateData = function(counts = NULL,isbatch = F,batchstr = 1,groupstr=1,numGenes=1000,propGrp=0.2,propBatch=0.5,gcoeffs=NULL,bcoeffs=NULL,myseed=10000){
  # counts is the matrix for reads count, half of treatment and half of control
  # batchstr is the strength of batch effect
  # groupstr is the strength of group effect
  suppressPackageStartupMessages(library(polyester));
  params = get_params(counts)
  batch = rep(c(1,-1),dim(counts)[2]/2)
  group = rep(c(1,-1),each=dim(counts)[2]/2)
  
  numSigGene = floor(numGenes * propGrp);
  
  if(is.null(gcoeffs)){
    gcoeffs = c(rnorm(numSigGene,sd=groupstr),rep(0,numGenes-numSigGene));  
  }
  
  
  numBatchGene = floor(numGenes * propBatch);
  tmp = c(rnorm(numBatchGene,sd=batchstr),rep(0,numGenes-numBatchGene));
  if(is.null(bcoeffs)){
    bcoeffs = bcoeffs = sample(tmp,numGenes,replace=F);
  }
  
  if(isbatch==F){
    bcoeffs=0*bcoeffs;
  }
  
  coeffs = cbind(bcoeffs,gcoeffs)
  controls = (bcoeffs != 0) & (gcoeffs==0)
  
  mod = model.matrix(~-1 + batch + group)
  dat0 = create_read_numbers(params$mu,params$fit,
                             params$p0,m=dim(counts)[1],n=dim(counts)[2],
                             beta=coeffs,mod=mod,seed=myseed);
  colnames(dat0) = colnames(counts);
  rownames(dat0) = 1:dim(dat0)[1];
  rm(mod);
  mod1 = model.matrix(~group)
  mod0 = cbind(mod1[,1])
  rownames(mod1) = colnames(counts);
  rownames(mod0) = colnames(counts);
  betamatrix = cbind(gcoeffs,bcoeffs)
  result = list("data"=dat0,"mod1"=mod1,"mod0"=mod0,"betamatrix"=betamatrix);
  return(result);
}

checkPCA = function(mydata=NULL,mincov = 10){
  tmpdata = mydata[rowSums(mydata)>mincov,]
  tmpdata = apply(tmpdata,2,function(x){x=x/sum(x)*1e6;return(x)})
  t_tmpdata = log2(t(tmpdata)+1);
  pca_proc <- prcomp(t_tmpdata[,apply(t_tmpdata, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE,retX=TRUE)
  load <- pca_proc$x[,1:2] 
  plot(load,pch=20,xlim=range(load[,1])*1.1,ylim=range(load[,2])*1.1) # set up plot 
  text(load*1.05,labels=colnames(mydata),cex=.7) # add variable names
  #biplot(pca_proc,xlim=c(-1,1),ylim=c(-1,1))
}


plotCoverage = function(x,treatmentCols=NULL,controlCols=NULL,genome=NULL,getFasta=F,cds.sort=cds.sort){
  
  # coordinates: genomic position
  #x = sigsites[sigsites[,"region"] == "chr:1180678-1180834",];
  #cat(paste(x[1,],collapse="\n"));
#   if(x[1,"region"]=="chr:932351-932351"){
#     browser();
#   }
  x1 = as.integer(as.numeric(x[,"start.b"]));
  minx1 = min(x1);
  maxx1 = max(x1);
  regionsize = maxx1-minx1+1;
  
  if(regionsize<200){
    regionsize=200;
    maxx1 = as.integer(maxx1+regionsize/2);
    minx1 = as.integer(minx1-regionsize/2);
  }
  
  x2 = minx1:maxx1;
  
  y = rep(0,length(x2)); # for control
  names(y)=as.character(x2);
  y[as.character(x1)]=rowMeans(x[,controlCols]);
  y[y<=0]=0;
  
  y1 = rep(0,length(x2)); # for treatment
  names(y1)=as.character(x2);
  y1[as.character(x1)]=rowMeans(x[,treatmentCols]);
  y1[y1<=0]=0;
  
  
  # significant sites
  sigsite2 = as.integer(as.numeric(x[x[,"score"]>=2,"start.b"]));
  sigsite3 = as.integer(as.numeric(x[x[,"score"]>=3,"start.b"]));
  sigsite5 = as.integer(as.numeric(x[x[,"score"]>=5,"start.b"]));
  
  # CDS sites
  if(getFasta==T & file.exists(genome)){
    #browser();
    #if(x[1,"region"]=="Ca21chr1:506609−507003"){
      
    #}
    
    myregion = data.frame("chr"=gsub("(.*):(.*)-(.*)","\\1",x[1,"region"]),"start"=minx1-1,"end"=maxx1,"name"=x[1,"region"],"score"=x[1,"score"],"strand"=x[1,"strand"],stringsAsFactors = F);
    
    cdssites=NULL;
    if(!is.null(cds.sort)){
      allsites = data.frame("chr"=gsub("(.*):(.*)-(.*)","\\1",x[1,"region"]),"start"=x2,"end"=x2+1,"name"=x[1,"region"],"score"=0,"strand"=x[1,"strand"],stringsAsFactors = F);
      
      tmpCDSsites = bedr(
        engine = "bedtools", 
        input = list(a = allsites  ,b=cds.sort[,1:6]), 
        method = "intersect", 
        params = "-s -loj -sorted",
        check.chr=F,
        verbose=F
      );
      cdssites = tmpCDSsites[tmpCDSsites[,"V4"]!=".","start"]
    }
    
    #myregion = data.frame("chr"="chr","start"=10000,"end"=10200,"name"="test","strand"="+",stringsAsFactors = F);
    #genome="~/gseq/prog/database/Genome/UCSC/B.subtilis/allchr.fa";
    tmp = get.fasta(myregion,fasta=genome,strand=T,check.chr=F,check.zero.based = F,verbose = F,check.merge = F);
    
    # get structure
    
    mycmd = paste("echo '",tmp[1,2],"' | RNAfold",sep="");
    mysst = system(mycmd,intern = T);
    mysst = mysst[2];
    mysst = gsub("(.*) (.*)","\\1",mysst);
    
    
    charAry = unlist(strsplit(tmp[1,"sequence"],split=""));
    charSSTAry = unlist(strsplit(mysst,split=""));
    
    charCol = rep("black",length(charAry));
    if(x[1,"strand"]=="-"){
      names(charCol) = as.character(x2)[length(x2):1];
    } else {
      names(charCol) = as.character(x2);
    }
    if(!is.null(cdssites)){charCol[as.character(cdssites)] = "darkgrey";};            
    charCol[as.character(sigsite2)] = trop[3];
    charCol[as.character(sigsite3)] = trop[4];
    charCol[as.character(sigsite5)] = trop[5]; 
    numCharPerLine = 60;
    numLine=floor(length(charAry)/numCharPerLine)+1;
    
    charAryToMatrix = function(charAry=NULL,numLine=NULL,numCharPerLine=NULL,defaultchar=""){
      newcharAry = rep(defaultchar,numLine*numCharPerLine);
      newcharAry[1:length(charAry)] = charAry;
      charMatrix = t(matrix(newcharAry,nrow=numCharPerLine));
      return(charMatrix);
    }
    charMatrix = charAryToMatrix(charAry = charAry,numLine=numLine,numCharPerLine=numCharPerLine);
    charSSTMatrix = charAryToMatrix(charAry = charSSTAry,numLine=numLine,numCharPerLine=numCharPerLine);
    charColMatrix = charAryToMatrix(charAry = charCol,numLine=numLine,numCharPerLine=numCharPerLine,defaultchar="white");
    #text(x = 0,y=numLine,paste(charAry[1:60],collapse=""),col="red",pos = 4)
    #textplot(0,numLine,paste(charAry[1:60],collapse=""),xlim=c(0,80))
    figwidth=20;
    figheight=numLine;
    
    pdf(paste(x[1,"region"],".fas.pdf",sep=""),width = figwidth,height=figheight);
    #plot(-1000,-1000,ylim=c(-1,numLine),xlim=c(-5,numCharPerLine+5),type="n", axes=FALSE, xlab="", ylab="",main=x[1,"region"]);
    plot(-1000,-1000,ylim=c(-1,numLine*2),xlim=c(-5,numCharPerLine+5),type="n", axes=FALSE, xlab="", ylab="",main=paste(myregion[1,"chr"],":",myregion[1,"start"],"-",myregion[1,"end"],":",myregion[1,"strand"],sep=""));
    for(i in 1:numLine){
      ypos = (numLine-i)*2;
      ypos1 = ypos - 1;
      mychars = charMatrix[i,];
      mySSTchars = charSSTMatrix[i,];
      mycols = charColMatrix[i,];
      newplot = ifelse(i==1,T,F);
      textplot(1:numCharPerLine,rep(ypos,numCharPerLine),words=mychars,cex = 1.5,col=mycols,new=F);
      textplot(1:numCharPerLine,rep(ypos1,numCharPerLine),words=mySSTchars,cex = 1.2,col=mycols,new=F);
    }
    
    dev.off();
    
  }
  
  # plot
  options( digits = 3 );
  title = paste(x[1,"region"],"; Gene: ",x[1,"gene"]," (","Evalue: ",format(x[1,"evalue"]),"; ","WeightedFold: ",format(x[1,"weightFold"]),")",sep="");
  figwidth=16/200*length(x2);
  figheight=figwidth*0.5;
  pdf(paste(x[1,"region"],".pdf",sep=""),width = figwidth,height=figheight);
  plot(x2,y,type='h',ylim=c(range(c(y,y1))),col=adjustcolor(trop[2],alpha.f=0.5),lwd=4,main=title,ylab="Normalized Coverage(log2)",xlab="Chromosome Position");
  lines(x2,y1,type='h',ylim=c(range(c(y,y1))),col=adjustcolor(trop[1],alpha.f=0.5),lwd=4);
  #abline(v=x3,lty=2,lwd=0.5,col="grey");
  Arrows(sigsite2,rep(-0.05,length(sigsite2)),sigsite2,rep(0,length(sigsite2)),code = 2,col=trop[3],arr.adj = 1,lwd=1);
  Arrows(sigsite3,rep(-0.05,length(sigsite3)),sigsite3,rep(0,length(sigsite3)),code = 2,col=trop[4],arr.adj = 1,lwd=1);
  Arrows(sigsite5,rep(-0.05,length(sigsite5)),sigsite5,rep(0,length(sigsite5)),code = 2,col=trop[5],arr.adj = 1,lwd=1);
  legend("topright",legend=c("treatment","control","P<1e-2","P<1e-3","P<1e-5"),lty=rep(1,5),lwd=c(2,2,4,4,4),col=trop);
  grid(col="grey",lwd=1.5);
  dev.off();
  options( digits = 7 );
  
  # get fasta
  
  
  return(T);
}


annoBioType = function(sigregion=NULL,genebed = NULL){
  # this function is designed for oranism without intron, only one exon for each gene
  sigMidpoint = round((sigregion[,2] + sigregion[,3])/2,0);
  sigregion[,2] =  sigMidpoint;
  sigregion[,3] =  sigMidpoint+1;
  sigregion.sort = bedr.sort.region(x=sigregion,check.zero.based = F,check.chr = F,check.valid = F,check.merge = F,verbose = F);
  
  genebed.sort = bedr.sort.region(x=genebed,check.zero.based = F,check.chr = F,check.valid = F,check.merge = F,verbose = F);
  
  upregion = bedr(
    engine = "bedtools", 
    input = list(a = sigregion.sort[,1:6] ,b=genebed.sort[,1:6]), 
    method = "closest", 
    params = "-s -D b -d -io -id",
    check.chr=F,
    verbose=F
  );
  upregion[,"V13"] = as.numeric(upregion[,"V13"]);
  upregion = upregion[upregion[,"V10"]!="." & abs(upregion[,"V13"])<500,];
  regionanno=NULL;
  if(dim(upregion)[1]>0){
    regionanno = upregion[,c("V4","V10","V13")];
    colnames(regionanno) = c("region","regionType","closestCDSDistance");
    regionanno[,"regionType"] = paste("5UTR__",regionanno[,"regionType"],sep="");
  }

  innerregion = bedr(
    engine = "bedtools", 
    input = list(a = sigregion.sort[,1:6],b=cds.sort[,1:6]), 
    method = "intersect", 
    params = "-s -loj -sorted",
    check.chr=F,
    verbose=F
  );
  
  if(!is.null(regionanno)){
    innerregion = innerregion[innerregion[,"V4"]!="." & is.na(match(innerregion[,"region"],regionanno[,"region"])),];  
  } else {
    innerregion = innerregion[innerregion[,"V4"]!=".",];
  }
    
  if(dim(innerregion)[1]!=0){
    innerregion[,c("V2","V3","V5","start","end")] = apply(innerregion[,c("V2","V3","V5","start","end")],2,as.numeric);
    tmpdistance = ifelse(innerregion[,"V6"]=="+",innerregion[,"start"]-innerregion[,"V2"]+innerregion[,"V5"],innerregion[,"V3"]-innerregion[,"end"]+innerregion[,"V5"]);
    tmp = data.frame("region" = innerregion[,"region"],"regionType"=innerregion[,"V4"],"closestCDSDistance"=tmpdistance,stringsAsFactors = F);
    tmp[,"regionType"] = paste("CDS__",tmp[,"regionType"],sep="");
    regionanno = rbind(regionanno,tmp);
  }
    
  downregion = bedr(
    engine = "bedtools", 
    input = list(a = sigregion.sort[,1:6] ,b=last.cds.sort[,1:6]), 
    method = "closest", 
    params = "-s -D b -d -io -iu",
    check.chr=F,
    verbose=F
  );
  
  downregion[,"V13"] = as.numeric(downregion[,"V13"]);
  if(!is.null(regionanno)){
    downregion = downregion[downregion[,"V10"]!="." & abs(downregion[,"V13"])<500 & is.na(match(downregion[,"V4"],regionanno[,"region"])),];  
  } else {
    downregion = downregion[downregion[,"V10"]!=".",];
  }
    
  if(dim(downregion)[1]>0){
    tmp = downregion[,c("V4","V10","V13")];
    colnames(tmp) = c("region","regionType","closestCDSDistance");
    tmp[,"regionType"] = paste("3UTR__",tmp[,"regionType"],sep="");
    regionanno = rbind(regionanno,tmp);
  }
    
  if(is.null(regionanno)){
    regionanno = data.frame("region"="","regionType"="","closestCDSDistance"=0);
  }
  
  sigregion.bed.merge = merge(regionanno,sigregion,by="region",all.y=T);
  sigregion.bed.merge[,"cdsgene"] = gsub("(.*)__(.*)","\\2",sigregion.bed.merge[,"regionType"]);
  sigregion.bed.merge = merge(sigregion.bed.merge,cds.length,by="cdsgene",all.x=T);
  sigregion.bed = sigregion.bed.merge[,c("chr","start","end","region","evalue","strand","gene","weightFold","regionType","closestCDSDistance","CDSLength")]
  return(sigregion.bed);
}


findClosetGene = function(sigregion=NULL,genebed = NULL){
  # this function is designed for oranism without intron, only one exon for each gene
  sigMidpoint = round((sigregion[,2] + sigregion[,3])/2,0);
  sigregion[,2] =  sigMidpoint;
  sigregion[,3] =  sigMidpoint+1;
  sigregion.sort = bedr.sort.region(x=sigregion,check.zero.based = F,check.chr = F,check.valid = F,check.merge = F,verbose = F);
  genebed = genebed[!is.na(match(genebed[,1],sigregion[,1])),];
  genebed.sort = bedr.sort.region(x=genebed,check.zero.based = F,check.chr = F,check.valid = F,check.merge = F,verbose = F);
  
  upregion = bedr(
    engine = "bedtools", 
    input = list(a = sigregion.sort[,1:6] ,b=genebed.sort[,1:6]), 
    method = "closest", 
    params = "-s -D b -d ",
    check.chr=F,
    verbose=F
  );
  upregion[,"V13"] = as.numeric(upregion[,"V13"]);
  upregion = upregion[upregion[,"V10"]!="." & abs(upregion[,"V13"])<500,];
  regionanno=NULL;
  if(dim(upregion)[1]>0){
    regionanno = upregion[,c("V4","V10","V13")];
    colnames(regionanno) = c("region","gene","closestCDSDistance");
  }
  return(regionanno);
}
