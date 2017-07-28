READS=$1
CPUNUM=$2
REF=$3
LOG=$4
OUTPUT=$5
HALFCPU=$((CPUNUM/2+0))
zcat $READS | bowtie2 -p $CPUNUM --no-unal -k 30  --fast -x $REF - 2> $LOG | CPUNUM=${HALFCPU} ~/prog/script/parse_bam_best_parallel_random.sh | samtools view -@ ${HALFCPU} -S -b -  > $OUTPUT
