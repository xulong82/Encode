# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - define the functions
# Rev: March 13, 2014

#---UCSC genomic features ---------------------------------------------------------------
setwd("~/Dropbox/Encode/UCSC")
chrom.info <- read.delim("./chromInfo.txt", header = F)  # chromosomes length
interval <- 200  # the binary genome unit size: the nucleotide length
chrom.length <- chrom.info$V2 %/% interval + 1  # chromosomes length in 200bp unit
genome.seg <- rep(0, sum(chrom.length))

#----------------------------------------------------------------------------------------
mySeg <- function(x) {  # Segment ENCODE broadPeak file
  # x: ENCODE broadPeak file
  # y: Segmentized broadPeak file with 200bp unit
  y <- genome.seg  # Initialize the return genome
  x$chrom <- gsub("chrX", "chr20", x$chrom)
  x$chrom <- gsub("chrY", "chr21", x$chrom)  # Encode has no Y chromosome
  x$chrom <- as.numeric(gsub("chr", "", x$chrom))
  pb <- txtProgressBar(min = 0, max = 100, style = 3)  # Progress bar  
  for (i in 1:nrow(x)) {
    if (i %% 1000 == 0) {
      progress <- (i * 100) %/% nrow(x)
      setTxtProgressBar(pb, progress)
    }  # Update progress bar
    x1 <- x[i, ]  # Single record
    offset <- ifelse(x1$chrom == 1, 0, sum(chrom.length[1:(x1$chrom - 1)]))
    start <- offset + x1$chromStart %/% interval
    end <- offset + x1$chromEnd %/% interval
    y[start:end] <- x1$V7
  }
  close(pb)  # Close progress bar
  return(y)  # genome wide binary signal of the input feature
}
#----------------------------------------------------------------------------------------