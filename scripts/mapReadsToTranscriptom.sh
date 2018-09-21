#!/usr/bin/env bash

READS=$1
CPUNUM=$2
REF=$3
LOG=$4
OUTPUT=$5
TMPID=$6
HALFCPU=$((CPUNUM/2+0))
TIRDCPU=$((CPUNUM/3+0))
QUARCPU=$((CPUNUM/4+0))
FIFTHCPU=$((CPUNUM/5+0))
echo  "zcat $READS | bowtie2 -p $HALFCPU --no-unal -k 30  --local -x $REF - 2> $LOG | CPUNUM=${FIFTHCPU} parse_bam_best_parallel_random.sh | samtools view -h -@ ${FIFTHCPU} -S -b -  | samtools sort -@ ${FIFTHCPU} -T ${TMPID} -m 1G -O BAM - > $OUTPUT"
eval "zcat $READS | bowtie2 -p $HALFCPU --no-unal -k 30  --local -x $REF - 2> $LOG | CPUNUM=${FIFTHCPU} parse_bam_best_parallel_random.sh | samtools view -h -@ ${FIFTHCPU} -S -b -  | samtools sort -@ ${FIFTHCPU} -T ${TMPID} -m 1G -O BAM - > $OUTPUT"
