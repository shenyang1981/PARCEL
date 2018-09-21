#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"));
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("bedr"))

parser <- ArgumentParser(description='extract the coverage information for significant regions');

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-o", "--outdir", type="character", nargs=1, help="prefix of directory for all output", metavar="coveragePlot",required=T);

parser$add_argument("-i", "--infile", type="character", nargs=1, help="significant regions as Input", metavar="combined_nag100_output2_wfilters.txt",required=T);

parser$add_argument("-s", "--samplefile", type="character", nargs=1, help="sample information", metavar="sampleList_Bsub_buffered.txt",required=T);

parser$add_argument("-b", "--usedBatch", type="character", nargs =1, help="batch ID for comparison (defined in comparison batch column) samplelist file", metavar="control",required=T);

parser$add_argument("-t", "--treatment", type="character", nargs =1, help="name of treatment condition in samplelist file", metavar="sam",required=T);

parser$add_argument("-c", "--control", type="character", nargs =1, help="name of control condition in samplelist file", metavar="control",required=T);


parser$add_argument("--covfile", type="character", nargs=1, help="coverage data stored as variable 'v1all' in Rdata file", metavar="combined_v1all.Rdata",required=T);

parser$add_argument("--vsall", default="T", type="character", metavar="T",help="whether use other samples in sample list as control [default %(default)s]");

parser$add_argument("--ismerge", default="F", type="character", metavar="F",help="whether merge libs with same condition and replicates [default %(default)s]");

parser$add_argument("--extWinSize", default=60, type="double", metavar=20,help="window size for calculating the wFC[default %(default)s]");

parser$add_argument("--genome", default="F", type="character", metavar="transcriptome.fas",help="transcriptome in fasta file[default %(default)s]");


# get parameters from command line -----------------------------------


args <- parser$parse_args()

infile = args$infile;
outdir = args$outdir;
covfile = args$covfile;
samplefile = args$samplefile;
usedBatch = args$usedBatch;
cond2 = args$treatment;
cond1 = args$control;
flagVsall=as.logical(args$vsall);
ismerge = as.logical(args$ismerge);
ext.win.size = args$extWinSize;
replot = as.logical(args$replot);
getFasta = as.logical(args$getFasta);
genome = args$genome;

# load sample information -------------------------------------------------

sampleInfo = read.table(samplefile,header=T,sep="\t",stringsAsFactors = F);

sampleInfo = sampleInfo[sampleInfo[,"ComparisonBatch"] == usedBatch,];

if(ismerge){
  sampleInfo[,"LibID"] = paste(sampleInfo[,"Condition"],sampleInfo[,"Replicates"],sep="__");
}
tmpSampleInfo = unique(sampleInfo[,c("Condition","LibID")]);

if(flagVsall==T){
  tmpSampleInfo[tmpSampleInfo[,"Condition"]!=cond2,"Condition"] = cond1;  
}

treatmentCols = tmpSampleInfo[tmpSampleInfo[,"Condition"]==cond2,"LibID"];
controlCols = tmpSampleInfo[tmpSampleInfo[,"Condition"]==cond1,"LibID"];


# read candidate region and coverage data ---------------------------------

output2 = fread(infile,header=T,sep="\t");
if(dim(output2)[1]>0){
  load(covfile);
  v1all  = v1all[chr %in% output2[,chr],];
  alllibs = names(v1all);
  alllibs = alllibs[is.na(match(alllibs,c("chr","pos")))];
  otherCols = alllibs[is.na(match(alllibs,c(treatmentCols,controlCols)))];
  
  covcols = c("chr","start","end","geneID","E_value","maxPos","maxFC","winSize");
  regionInfo = output2[,covcols,with=F];
  regionInfo[,region:=paste(chr,":",start,"-",end,sep="")]
  
  
  # filter out low complexity region ---------------------
  regionbed = regionInfo[,c("chr","start","end","region","E_value"),with=F];
  tmpMidpoint = regionbed[,round((start + end)/2,0)];
  tmp.ext.win.size = 30;
  regionbed[,strand:= ifelse(gsub("(.*):(.*)","\\2",chr)=="Neg","-","+")];
  regionbed[,chr:= gsub("(.*):(.*)","\\1",chr)];
  regionbed[,end:= ifelse(strand=="+", tmpMidpoint + tmp.ext.win.size,tmpMidpoint)];
  regionbed[,start:= ifelse(strand=="-", tmpMidpoint - tmp.ext.win.size,tmpMidpoint)];
  regionbed = data.frame(unique(regionbed));
  regionbed.sort = bedr.sort.region(x = regionbed,engine = "bedtools",check.zero.based = F,check.chr = F,check.valid = F,check.merge = F);
  regionnucinfo = data.table(get.fasta(regionbed.sort,fasta = genome,strand = T,use.name.field = T,check.chr = F,check.valid = F,check.zero.based = F,check.sort = F,check.merge = F));
  regionnucinfo[,Lcomplex:= unlist(lapply(sequence,function(x){
    nucCount = table(unlist(strsplit(toupper(x),split="")));
    if(max(nucCount/sum(nucCount))>0.5){
      return(T);
    } else {
      return(F);
    }
  }))];
  regionInfo = regionInfo[!region%in%regionnucinfo[Lcomplex==T,index],];
  
  if(dim(regionInfo)[1]>0){
    
    # write filtered regions
    write.table(output2[region%in%regionInfo[,region],],file=paste(outdir,"combined_",cond2,"_output2_wfilters.noLComplex.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
    
    # extend candidate region on both side ---------------------------
    
    regionInfo[,c("extstart","extend"):=.(start-ext.win.size,end+ext.win.size)];
    regionInfo[,extstart:=ifelse(extstart>=1,extstart,1)];
    
    # coverage information were stored in v1all ----------------------
    
    v1all[,c("start","end"):=.(pos,pos+1)];
    setkeyv(v1all,c("chr","start","end"));
    setkeyv(regionInfo,c("chr","extstart","extend"));
    
    # find the coverage information using foverlaps ------------------
    
    covoutput = foverlaps(v1all,regionInfo[,c("chr","extstart","extend","region"),with=F],by.x=c("chr","start","end"),by.y=c("chr","extstart","extend"),nomatch = 0,which = F,type = "any");
    
    # fill zero to position without coverage information ---------------------- 
    
    fillzero = function(x){
      tmpcols = names(x)[c(-1,-2,-3)]; # pos, start and end
      tmprange = as.matrix(x[1,2:3])[1,] # get the range of position
      tmppos = data.table("pos" = seq(tmprange[1],tmprange[2],1), "extstart" = tmprange[1], "extend" = tmprange[2]);
      tmpcov = merge(tmppos,x,by="pos",all.x=T);
      tmpcov[,c(tmpcols) := lapply(.SD,function(x){x[is.na(x)]=0;return(x)}),.SDcols=tmpcols];
      return(tmpcov);
    }
    covoutput = covoutput[,fillzero(.SD),by=region,.SDcols=c("pos","extstart","extend",treatmentCols,controlCols,otherCols)];
    covoutput[,c("nucstart","nucend"):=.(pos,pos+1)];
    covoutput = merge(regionInfo,covoutput,by="region",all.y=T);
    covoutput[,relpos:=ifelse(grepl(":Pos",chr),pos-start+1,end-pos+1)];
    
    # extract nucleotide information from transcriptome fasta file used in mapping ----------
    
    candbed = as.data.frame(covoutput[,c("chr","nucstart","nucend"),with=F]);
    candbed[,"name"] = "tmp";
    candbed[,"score"] = 0;
    candbed[,"strand"] = ifelse(gsub("(.*):(.*)","\\2",candbed[,"chr"])=="Neg","-","+");
    candbed[,"chr"] = gsub("(.*):(.*)","\\1",candbed[,"chr"]);
    candbed = unique(candbed);
    candbed.sort = bedr.sort.region(x = candbed,engine = "bedtools",check.zero.based = F,check.chr = F,check.valid = F,check.merge = F);
    nucinfo = data.table(get.fasta(candbed.sort,fasta = genome,strand = T,use.name.field = F,check.chr = F,check.valid = F,check.zero.based = F,check.sort = F,check.merge = F));
    
    # merge covoutput with nucleotide 
    covoutput[, cleanchr := gsub("(.*):.*", "\\1", chr)];
    covoutput[,nucindex:=paste(cleanchr,":",nucstart-1,"-",nucend-1,"(",ifelse(gsub("(.*):(.*)","\\2",chr)=="Pos","+","-"),")",sep="")];
    covoutputWithNuc = merge(covoutput,nucinfo,by.x="nucindex",by.y="index");
    # preparing final output
    covoutputcols = c("geneID","region","chr","start","end","E_value","maxPos","maxFC","winSize","extstart","extend","nucstart","relpos","sequence",treatmentCols,controlCols,otherCols);
    covoutputWithNuc = covoutputWithNuc[,covoutputcols,with=F];
    covoutputWithNuc = unique(covoutputWithNuc[order(geneID,region,relpos),]);
    setnames(covoutputWithNuc,"sequence","Nucleotide");
    setnames(covoutputWithNuc,"nucstart","position");
    write.table(covoutputWithNuc,file=paste(outdir,"combined_",cond2,"_covinfo.xls",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
    
    # concatnate nucleotides to sequence
    seqinfo = covoutputWithNuc[,paste(Nucleotide,collapse = ""),by=region];
    names(seqinfo) = c("region","sequence");
    write.table(seqinfo,file=paste(outdir,"combined_",cond2,"_seqinfo.xls",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
  }
} else {
  covoutputcols = c("geneID","chr","start","end","E_value","maxPos","maxFC","winSize","extstart","extend","position","relpos","Nucleotide");
  nullresult = data.frame(matrix(ncol=length(covoutputcols)));
  colnames(nullresult) = covoutputcols;
  write.table(nullresult[!is.na(nullresult[,1]),],file=paste(outdir,"combined_",cond2,"_covinfo.xls",sep=""),col.names=T,row.names=F,sep="\t",quote=F);
}

