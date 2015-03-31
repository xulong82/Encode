# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - define the functions
# Rev: March 13, 2014

#---UCSC genomic features ---------------------------------------------------------------
setwd("~/Dropbox/Info")
chrom.info <- read.delim("./chromInfo.txt", header = F)  # chromosomes length
interval <- 200  # the binary genome unit size: the nucleotide length
chrom.length <- chrom.info$V2 %/% interval + 1  # chromosomes length in 200bp unit
genome.seg <- rep(0, sum(chrom.length))

refgene <- read.delim("myRefGene.tx", header = T)  # RefSeq transcripts
refexon <- read.delim("myRefGene.exons", header = F)  # RefSeq exons
refgene.pm <- data.frame(chrom      = refgene$chrom,  # Promoter: -2k -> +2k TSS 
                         chromStart = refgene$txStart - 2000,
                         chromEnd   = refgene$txStart + 2000)
refgene.tx <- data.frame(chrom      = refgene$chrom,  # Transcript
                         chromStart = refgene$txStart,
                         chromEnd   = refgene$txEnd)
refgene.ex <- data.frame(chrom    = refexon$V1,  # Exon
                         chromStart = refexon$V2,
                         chromEnd   = refexon$V3)

#----------------------------------------------------------------------------------------
