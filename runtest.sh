#!/bin/bash

# Run one simple test to check program installation

./segmentanalysis.py -f -v -s 30 examples/dm6.Xand2L.fa.gz :X:11982050:12772070 examples/DmelMapTable.Xand2L.bed
#./genomecorrelation.py -m -c X,2L examples/dm6.nounmapped.fa.gz.X-11982050-12772070/cytocounts.l30-11982050-merged.txt examples/Berlin.ectopics.tsv
