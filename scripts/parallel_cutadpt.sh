#!/usr/bin/env bash

INPUT=$1
OUTPUT=$2
APT1=$3
APT2=$4
THREADS=$5
OUTLOG=$6
GZFQ=${GZFQ:="yes"}
POLYA=${POLYA:="no"}
MINOVER=${MINOVER:=5}
QUALTRIM=${QUALTRIM:=""}
if [ "$GZFQ" == "yes" ] 

then
	PIGZ=" | pigz -p $THREADS -c "
else
	PIGZ=""
fi

if [ "$POLYA" == "yes" ]

then
	TRIMPOLYA=" -a POLYA=A\{100\} "
else
	TRIMPOLYA=""
fi


if [ "$QUALTRIM" != "" ]

then
	QUALTRIM=" -q $QUALTRIM "
else
	QUALTRIM=""
fi

echo "zcat $INPUT | parallel -j $THREADS --pipe -N1000000 --block 100M \"cutadapt ${QUALTRIM} -a AP1=${APT1} -g AP2=${APT2} ${TRIMPOLYA}  -n 2 -O ${MINOVER} -y 'FOUND__{name}__' -m 18 - \" $PIGZ > $OUTPUT 2>$OUTLOG"
eval "zcat $INPUT | parallel -j $THREADS --pipe -N1000000 --block 100M \"cutadapt ${QUALTRIM} -a AP1=${APT1} -g AP2=${APT2} ${TRIMPOLYA}  -n 2 -O ${MINOVER} -y 'FOUND__{name}__' -m 18 - \" $PIGZ > $OUTPUT" 2>$OUTLOG
