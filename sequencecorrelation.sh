#!/bin/bash
# Run cytocouns collection and correlation search for one sequence
set -e

fragment=${1}
./sequencecorrelation.py -v -e X:11982050:12772070 examples/dm6.nounmapped.fa.gz ${fragment} examples/DmelMapTable.160615c.bed
for corr in Pearson Spearman; do
    ./genomecorrelation.py -m -c X -r $corr  examples/dm6.nounmapped.fa.gz.${fragment}.cytocounts.txt examples/Berlin.ectopics.tsv -s ${fragment}.corr.png
done



