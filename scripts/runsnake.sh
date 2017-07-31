#!/usr/bin/env bash

SMAKE=$1
SCONF=$2
JOBNAME=$3
RUNTIME=$4
NUMJOBS=$5
CONFIGVAR=$6
if [ -z "$RUNTIME" ]; then
	RUNTIME="48"
fi
if [ -z "$NUMJOBS" ]; then
	NUMJOBS="8"
fi
if [ -z "$CONFIGVAR" ]; then
	CONFIGVAR='DUMPVAR=NULL'
fi
# snakemake requires python3
source activate snakemake-tutorial;
# qsub for snakemake itself
qsub="qsub -pe OpenMP 1 -l mem_free=1G -l h_rt=48:00:00 -m bes -j y -V -b y -cwd";
# -j in cluster mode is the maximum number of spawned jobs

# here doc only for debugging. otherwise py3k not available
#cat <<EOF
$qsub -N ${JOBNAME}.snakemake -o ${JOBNAME}.snakemake.qsub.log "snakemake --jn \"${JOBNAME}.skj.{rulename}.{jobid}.sh\" -j ${NUMJOBS} -c \"qsub -pe OpenMP {threads} -l mem_free=50G -l h_rt=${RUNTIME}:00:00 -j y -V -b y -cwd\" -s ${SMAKE} --configfile ${SCONF} --config ${CONFIGVAR} -w 60";
#EOF
