#!/bin/sh
#
# RNA-seq data analysis pipeline 
# Author: XuLong Wang (xulong.wang@jax.org)

# Read count estimation with htseq-count

echo $0
echo "begin: `date`"

dir_sou="/data/xwang/Encode/BAM"
dir_des="/data/xwang/Encode/BAM"
mygtf="/data/xwang/GTF/mm9Gene.gtf"
#
FILES=`find $dir_sou -name *.bam`

for filename in $FILES; do
  echo $filename
  filename2=$(basename $filename)
  echo $filename2
#  samtools view -h "$dir_sou"/$filename2 > "$dir_des"/${filename2/.bam/.sam}
   htseq-count -s reverse "$dir_sou"/${filename2/.bam/.sam} $mygtf > "$dir_des"/${filename2/.bam/.txt}
done

echo "end: `date`"

