
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("reshape2"));
suppressPackageStartupMessages(library("data.table"));
suppressPackageStartupMessages(library("Matrix"));
parser <- ArgumentParser(description='conversion between "wide" and "long"');

# specify our desired options
# by default ArgumentParser will add an help option
#parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
#parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")

parser$add_argument("-i", "--infile", type="character", nargs=1, help="txt (zipped) file as Input", metavar="mytable.txt.gz",required=T);

parser$add_argument("-o", "--outfile", type="character", nargs=1, help="txt (zipped) file as Output", metavar="mytable_convert.txt.gz",required=T);

parser$add_argument("-m", "--measure", type="character", default="", help="measure.vars, where each element of the list contains the columns that should be combined together. Group of variable is separated by :", metavar="var1,var2,var3:var4,var5,var6");

parser$add_argument("--header", default=T, type="character", metavar="T",help="Whether there is header in input file [default %(default)s]");

parser$add_argument("--natozero", default=T, type="character", metavar="T",help="NA would be replaced with 0 [default %(default)s]");

parser$add_argument("--filterbySum", default=0, type="double", metavar="0",help="select rows by sum of row [default %(default)s]");

parser$add_argument("--downSampleProp", default=1, type="double", metavar="1",help="downsample coverage reads [default %(default)s]");

parser$add_argument("--pigzthreads", default=6, type="double", metavar="6",help="downsample coverage reads [default %(default)s]");

parser$add_argument("--sampleinfo", default="", type="character", metavar="sampleinfo.txt",help="merge columns from same condition[default %(default)s]");

#parser$add_argument("--batch", default="", type="character", default="",metavar="batch1",help="which comparison batch should be used");

parser$add_argument("-v", "--variable", type="character", default="", help="variable name for merged variables", metavar="metavar1:metavar2",required=T);

parser$add_argument("-s", "--splitBy", type="character", default="", help="variable name for splitting datasets", metavar="id");

parser$add_argument("-n", "--names", type="character", default="", help="variable names for label variables", metavar="V1,V2");

parser$add_argument("-f", "--formula", type="character", default="", help="formula for dcasting long format", metavar="id1+id2~var3");

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

# print some progress messages to stderr if "quietly" wasn't requested
#if ( args$verbose ) {
#  write("writing some verbose output to standard error...\n", stderr())
#}

#setwd("~/gseq/prog//Chengqi/metaAnalysis/")

iname = args$infile;
output = args$outfile;
measureVars = args$measure;
valueVars = args$variable;
labelVars = args$names;
myformula = args$formula;
isheader = as.logical(args$header);
replaceNA = as.logical(args$natozero);
filterbySum = args$filterbySum;
sampleinfo = args$sampleinfo;
downSampleProp = as.numeric(args$downSampleProp);
splitBy = args$splitBy;
pigzthreads = args$pigzthreads;
#batch = args$batch;


pigz_pipe = function(filename, mode="read", cores=4) {
  if(mode == "read") {
    con <- pipe(paste0("cat ", filename, " | pigz -dcp ", cores), "rb")
  } else {
    con <- pipe(paste0("pigz -p ", cores, " > ", filename), "wb")
  }
  con
}

downSample = function(x=NULL,prop = NULL){
  y = rbinom(n = length(x),size=x,p=prop);
  return(y);
}

cat(paste(c(iname,output,measureVars,valueVars,myformula,isheader,"\n"),collapse=" ; "));

#iname = "melt_enhanced.csv"
# iname = "test.txt.gz"
# output = "testout.txt.gz"
# isheader=F;
# valueVars = "V3";
# myformula = "V1+V2~V4";

if(grepl("gz$",iname)){
  #iname = paste("zcat ",iname,sep="");
  iname = paste("gzip -cd ",iname,sep="");
}
DT = fread(iname,header = as.logical(isheader));

processDT = function(DT=NULL,valueVars=NULL,measureVars=NULL,myformula=NULL,replaceNA=NULL,filterbySum=NULL,labelVars=NULL){
  # the labelVars are used while measureVars and myformula are not provided
  myvars = unlist(strsplit(valueVars,split=','));
  if(measureVars!=""  && myformula==""){
    cat("melt DT\n");
    #measureVars = lapply(unlist(strsplit("dob_child1,dob_child2,dob_child3:gender_child1,gender_child2,gender_child3",split=":")),function(x){unlist(strsplit(x,split=','))});
    measureVars = lapply(unlist(strsplit(measureVars,split=":")),function(x){unlist(strsplit(x,split=','))});
    outDT = melt(DT, measure = measureVars, value.name = myvars);  
  } else if (measureVars==""  && myformula!="") {
    cat("dcast DT\n");
    outDT = dcast(DT,myformula , value.var = myvars,fun.aggregate=sum);
  } else if (measureVars=="" && myformula=="") {
    cat("nothing to be done with DT\n");
    outDT = DT;
  } else {
    stop(sprintf("--measure and --formula can not be specified simultaneously!\n"));
  }
  #outDT = as.data.frame(outDT);
  
  # replace NA to 0
  
  if(replaceNA){
    cat("Replace NA\n");
    # replaceNA = function(DT=NULL) {
    #   for (i in names(DT))
    #     DT[is.na(get(i)),(c(i)):=0,with=FALSE]
    #   return(DT);
    # }
    #outDT=replaceNA(outDT);
    outDT = outDT[,lapply(.SD,function(x){x[is.na(x)]=0;return(x)})];
  }
  
  if(downSampleProp!=1){
    cat("Downsample DT\n");
    tmpvar = gsub(".*~(.*)","\\1",myformula);
    allcols = as.character(as.matrix(unique(DT[,tmpvar,with=F])));
    #print(allcols);
    for(tmpcol in allcols){
      outDT[[tmpcol]] = downSample(x = outDT[[tmpcol]],prop=downSampleProp);
    }
  }
  
  if(filterbySum!=0){
    cat("Filter DT by rowSum\n");
    if(myformula != ""){
      tmpvar = gsub(".*~(.*)","\\1",myformula);
      if(length(tmpvar)==1){
        #print(tmpvar);
        print(unique(DT[,tmpvar,with=F]));
        outDT=outDT[outDT[,Reduce('+',.SD),.SDcols=unique(DT[[tmpvar]])]>=filterbySum];  
      }  
    }
    if(labelVars != ""){
      ids = unlist(strsplit(labelVars,split=','));
      allids = names(outDT);
      valueids = allids[is.na(match(allids,ids))];
      outDT = outDT[outDT[,Reduce('+',.SD),.SDcols = valueids]>=filterbySum,];
    }
  }
  
  return(outDT);
}

# do it by each subset

if(splitBy!="" & splitBy!="NULL"){
  allvars = as.character(as.matrix(unique(DT[,splitBy,with=F])));  
  outDT=NULL;
  for(tmpvar in allvars){
    print(paste("Processing ",tmpvar,sep=":"));
    subDT = DT[DT[[splitBy]]==tmpvar,];
    subOutDT = processDT(DT=subDT,valueVars = valueVars,measureVars = measureVars,myformula=myformula,replaceNA=replaceNA,filterbySum=filterbySum,labelVars = labelVars);
    if(is.null(outDT)){
      outDT=subOutDT;
    } else {
      outDT=rbindlist(list(outDT,subOutDT), use.names=TRUE,fill=T);
    }
  }
} else {
  outDT = processDT(DT=DT,valueVars = valueVars,measureVars = measureVars,myformula=myformula,replaceNA=replaceNA,filterbySum=filterbySum,labelVars = labelVars);
}


# downsampling


# merge cols using sampleinfo

if(sampleinfo !="" & file.exists(sampleinfo)){
  sampleinfo = read.table(sampleinfo,header=T,sep="\t",stringsAsFactor=F);
  usedcols = colnames(outDT);
  sampleinfo = sampleinfo[!is.na(match(sampleinfo[,"LibID"],usedcols)),];
  conds = paste(sampleinfo[,"Condition"],sampleinfo[,"Replicates"],sep="__");
  allcols = sampleinfo[,"LibID"];
  for(mycon in unique(conds)){
    print(paste("Merging ",mycon,sep=":"));
    mycols = allcols[conds==mycon];
    outDT[, (c(mycon)) := rowSums(.SD,na.rm=T),.SDcols=mycols];
  }
}



print("output result\n");
if(grepl("gz$",output)){
  #gz1 <- gzfile(output, "w");
  pigzcon = pigz_pipe(output,mode="write",cores=pigzthreads);
  write.table(outDT,pigzcon,col.names=T,row.names=F,sep="\t",quote=F);
  close(pigzcon);
} else {
  write.table(outDT,file=output,col.names=T,row.names=F,sep="\t",quote=F);
}
