#!/bin/sh
#
# RNA-seq data analysis pipeline
# Author: XuLong Wang (xulong.wang@jax.org)

module load bowtie/1.0.0
module load rsem/1.2.16

echo $0
begin=`date +%h`
 
dir1="/hpcdata/xwang/Network/FASTQ2"
dir2="/hpcdata/xwang/Network/RSEM"
myref="/hpcdata/xwang/Rsem/GRCm38"

files=`find $dir1 -name '*.fastq'`

for name1 in $files; do
  name2=`basename $name1`
  name3=${name2/.fastq/}
  echo $name3
  rsem-calculate-expression -p 20 \
  			    --solexa-quals \
			    --forward-prob 0 \
                            "$dir1"/"$name3".fastq \
  			    "$myref" \
  			    "$dir2"/"$name3"
done

end=`date +%h`
echo $((end-begin))

# Options
# Make the genome bam (Caution: This take 10 hours per sample)
#        		    --output-genome-bam 
# Generating a Wiggle file
#rsem-bam2wig sorted_bam_input wig_output.wig wiggle_name

# Prepare reference sequences with RSEM
#rsem-prepare-reference --gtf ./GTF/mm10knownGene.gtf \
#		       --transcript-to-gene-map ./GTF/mm10knownIsoforms.txt \
#                       ./Genome/mm10.fa \
#                       ./RSEM/mm10rsem
#
#rsem-prepare-reference --gtf ./GTF/Mus_musculus.GRCm38.73.gtf \
#                       ./Genome/GRCm38.73 \
#                       ./RSEM/NCBIM37.59
#
#rsem-prepare-reference --gtf ./GTF/Mus_musculus.NCBIM37.59.gtf \
#                       ./Genome/mm9/regular \
#                       ./RSEM/NCBIM37
#
