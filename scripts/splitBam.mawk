#!/usr/bin/env bash

mawk -F "\t" -v blockSize=$BLOCKSIZE -v splitTag=$SPLITTAG -- 'BEGIN{
	OFS="\t";
	blockSize = blockSize>0?blockSize:500000;
	splitTag = splitTag!=""?splitTag:"--SplitTag--";
	lastid="";count=0}
{	
	if($0~ /^@/){print;} 
	else {
		currentid=$1;
		if( count>blockSize && lastid != currentid){
			print splitTag;
			count=0;		
		}
		print;
		count++;
		lastid=currentid;
	}
}'
