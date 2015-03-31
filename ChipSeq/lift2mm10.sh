#!/bin/bash

# Lift over ENCODE broadPeak from mm9 to mm10

echo $0

module load liftover

cd /home/xwang/Dropbox/Network/ChipSeq

dir1="/home/xwang/Dropbox/Network/ChipSeq/broadPeak"
dir2="/home/xwang/Dropbox/Network/ChipSeq/bed"
dir3="/home/xwang/Dropbox/Network/ChipSeq/mm10"

files=`find $dir1 -name *.broadPeak`

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/.broadPeak/.bed}
  name4=${name3/wgEncodeLicr/}
  name5=${name4/C57bl6StdPk/}
  echo $name5
  awk '{print $1"\t"$2"\t"$3"\t"$7}' $name1 > ${dir2}/${name5}
  liftOver ${dir2}/${name5} mm9ToMm10.over.chain ${dir3}/${name5} unmapped
done

