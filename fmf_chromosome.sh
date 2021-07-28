#!/bin/bash

set -e

if [ $# -ne 3 ]; then
    echo "This script runs segmentanalysis.py for particular chromosome of whole genome"
    echo "./fmf_chromosome.sh <genome> <cytobands> <chromosome>"
    exit
fi

genome=$1
bands=$2
chr=$3

sizes="30,50"
excludefr="repeats.txt"
includefr="include.txt"
threads=$(fgrep -c processor /proc/cpuinfo)

runstr="./segmentanalysis.py -f -s $sizes -b $excludefr -i $includefr $genome "

echo "Running parallel calculations of FMF frequencies"

echo "" > runlist.sh
cat $bands | awk -v FS="\t" -v OFS="\t" -v chr=$chr -v runstr="${runstr}" -v bands=$bands\
            '($1 == chr) {print runstr ":" chr ":" $2 ":" $3 " " bands}' >> runlist.sh

# Execute analysis in parallel
xargs --arg-file=runlist.sh --max-procs=$threads --replace /bin/sh -c "{}"

