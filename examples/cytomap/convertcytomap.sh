#!/bin/bash
cat DmelMapTable.160615c.csv | grep -v '#' | sed 's/:/\t/;s/\.\./\t/' | awk -v FS="\t" -v OFS="\t" '$3 != ""  {print $3,$4,$5,$3$1}' > DmelMapTable.160615c.bed
