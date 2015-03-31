#!/bin/sh
#
# RNA-seq data analysis pipeline 
# Author: XuLong Wang (xulong.wang@jax.org)

# Read trimming with Trimmomatic 

echo $0
begin=`date +%m`

module load java/1.7.0

mydir="/hpcdata/xwang/Encode/FASTQ2"
files=`find $mydir -name '*.fastq'`
#
for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/.fastq/}
  echo $name3
  java -jar /home/xwang/bin/Trimmomatic-0.32/trimmomatic-0.32.jar SE -threads 20 -phred64 \
            "$mydir"/"$name3".fastq \
            "$mydir"/"$name3"_cut.fastq \
            HEADCROP:6
done
#
#	    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
# Note: HEADCROP is how ENCODE does, LEADING...seems reasonable but gives strange results
#
end=`date +%m`
echo $((end-begin))
