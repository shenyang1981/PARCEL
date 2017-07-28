#!/usr/bin/env bash

if [ -z "$CPUNUM" ]; then
	CPUNUM="8"
fi

if [ -z "$BLOCKSIZE" ]; then
	BLOCKSIZE="500000"
fi
BLOCKSIZE=$BLOCKSIZE splitBam.mawk | SPLITTAG="--SplitTag--" parallel -j $CPUNUM --block 400M --recstart "--SplitTag--" --pipe -N1 "parse_bam_best_random.pl"
