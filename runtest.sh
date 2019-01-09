#!/bin/bash

# Run one simple test to check program installation

./segmentanalysis.py -v -s 30 -d 2.5 examples/dm6.nounmapped.fa.gz :X:11982050:12772070 examples/DmelMapTable.160615c.bed 
./genomecorrelation.py -m -c X,2L examples/dm6.nounmapped.fa.gz.X-11982050-12772070/cytocounts.l30-merged.txt examples/Berlin.ectopics.tsv 


