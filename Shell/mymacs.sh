#!/bin/sh
#
# ChIP-seq data analysis pipeline
# Author: XuLong Wang (xulong.wang@jax.org)

# Usage: sh mymacs.sh > mylog

module load python/2.7.3
module load MACS/1.4.2
module load bedtools/2.17.0

echo $0
echo "begin: `date`"

cd /data/xwang/Prdm9/Encode

for repID in 1 2; do
  macs14 -t wgEncodeLicrHistoneTestisH3k04me1MAdult8wksC57bl6StdRawDataRep"$repID".sam \
         -c wgEncodeLicrHistoneTestisInputMAdult8wksC57bl6StdRawDataRep"$repID".sam \
         -n EncodeTestisH3k04me1Rep"$repID" \
         -f SAM -g mm -s 36
done

# Intersecting the peaks of the replicates
intersectBed -a EncodeTestisH3k04me1Rep1_peaks.bed -b EncodeTestisH3k04me1Rep2_peaks.bed > EncodeTestisH3k04me1_peaks.bed

echo "end: `date`"

# echo "." | mail -s "Task in Rockhopper have completed!" emailofx@gmail.com 

