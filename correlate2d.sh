#!/bin/bash

#set -e

genome=$1
chr=$2
chrsize=$3

chunk=10000
threads=$(fgrep -c processor /proc/cpuinfo)

start=0
i=0
echo "" > runlist.sh

# Generate regions
while [ $start -lt $chrsize ]; do
    
      let "stop = start + chunk"
      if [ $stop -gt $chrsize ]; then
         let "stop = chrsize - 1 "
      fi
      
      ii=$( printf %06d $i )
      echo "./segmentanalysis.py -f --nodump -b blacklist-n.seq -s 50 $genome :$chr:$start:$stop && mv $genome.$chr-$start-$stop/rawcounts.l50-$start-merged.txt $genome.$chr.l50-$ii-merged.txt ; rm -rf $genome.$chr-$start-$stop" >> runlist.sh
      let "start = start + chunk"
      let "i = i + 1"
done

# Execute analysis in parallel
#xargs --arg-file=runlist.sh --max-procs=$threads --replace /bin/sh -c "{}"


start=0
i=0
echo -e "Start\tStop\tCount" > $genome.$chr.seqcorr.txt

while [ $start -lt $chrsize ]; do
      let "stop = start + chunk"
      if [ $stop -gt $chrsize ]; then
         let "stop = chrsize - 1 "
      fi
      
      ii=$( printf %06d $i )

      if [ -f $genome.$chr.l50-$ii-merged.txt ]; then
          awk -v stop=$stop -v chunk=$chunk -v OFS="\t" '{print stop, $2+chunk, $5}' $genome.$chr.l50-$ii-merged.txt  >> $genome.$chr.seqcorr.txt
      fi
      let "start = start + chunk"
      let "i = i + 1"

done
