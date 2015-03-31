#!/bin/sh
# Use ChromHMM to find chromatin states
# Author: Xulong Wang (xulong.wang@jax.org)

HMMDIR="/Users/xwang/Applications/ChromHMM"
DATADIR="/Users/xwang/Dropbox/Encode"

usage() { 
    echo "Usage: `basename $0` arg1 [1|2]" 
}

if [ $# -ne 1 ]; then
    usage
    exit 1
fi

EXE=$1

echo "begin: `date`"

if [ $EXE -eq 1 ]; then
    # Binarize bed files
    echo "Start binarizing files ..."
    java -mx1600M -jar $HMMDIR/ChromHMM.jar \
        BinarizeBed -b 200 -peaks \
            $HMMDIR/CHROMSIZES/mm10.txt \
            $DATADIR/broadPeak \
            $DATADIR/HMM/cellmarkfiletable \
            $DATADIR/HMM/binary
fi

if [ $EXE -eq 2 ]; then
    # Learn Model
    echo "Start learning models ..."
    java -mx1600M -jar $HMMDIR/ChromHMM.jar \
        LearnModel -b 200 \
            $DATADIR/HMM/binary \
            $DATADIR/HMM/model 5 mm10
fi

echo "end: `date`"

