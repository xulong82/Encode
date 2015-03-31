# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - define the functions
# Rev: March 13, 2014

library(DESeq)

#---UCSC genomic features ---------------------------------------------------------------
setwd("~/Dropbox/Encode/UCSC")
chrom.info <- read.delim("./chromInfo.txt", header = F)  # chromosomes length
interval <- 200  # the binary genome unit size: the nucleotide length
chrom.length <- chrom.info$V2 %/% interval + 1  # chromosomes length in 200bp unit
genome.binary <- rep(0, sum(chrom.length))

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
myCount2 <- function(x1, x2) {  # Calculate number of 1 in a region
  # x1: BED columns 1-3 (1:chrom; 2:chromStart; 3:chromEnd)
  # x2: Binary genome/epigenome features
  # y:  number of 1 in x2 with regions defined by x1
  y <- c(rep(0, length(x1)))
  x1$chrom <- gsub("chrX", "chr20", x1$chrom)
  x1$chrom <- gsub("chrY", "chr21", x1$chrom)  # Encode has no Y chromosome
  x1$chrom <- as.numeric(gsub("chr", "", x1$chrom))
  pb <- txtProgressBar(min = 0, max = 100, style = 3)  # Progress bar  
  for (i in 1:nrow(x1)) {
    if (i %% 1000 == 0) {
      progress <- (i * 100) %/% nrow(x1)
      setTxtProgressBar(pb, progress) 
    }  # Update progress bar
    x1.1 <- x1[i, ]  # Single record
    offset <- ifelse(x1.1$chrom == 1, 0, sum(chrom.length[1:(x1.1$chrom - 1)]))
    start <- offset + x1.1$chromStart %/% interval
    end <- offset + x1.1$chromEnd %/% interval
    y[i] <- sum(x2[start:end])
  }
  close(pb)  # Close progress bar
  return(y)
}
#----------------------------------------------------------------------------------------
myQnorm <- function(x) {  # Quantile normalization to normal distribution
  # x: y: 
  norm <- rnorm(length(x), mean = median(x), sd = sd(x) / sqrt(length(x)))
  y <- sort(norm)[rank(x)]
  return(y)
}
#----------------------------------------------------------------------------------------
myVst <- function(x1, x2) {
  # x1: data frame of read count, x2: conditions
  cds <- newCountDataSet(x1, x2)
  cds <- estimateSizeFactors(cds)  # Library size normalization
  cds <- estimateDispersions(cds)  # Dispersion estimation
  cdsblind <- estimateDispersions(cds, method = "blind")
  vsd <- getVarianceStabilizedData(cdsblind)
  return(vsd)
}
#----------------------------------------------------------------------------------------
# mySearch <- function(x1, x2) {  # Search overlaps between genomic features
#   # x1: broadPeak files; x2: ucsc.genes | ucsc.exons
#   y <- matrix(0, nrow = nrow(x1), ncol = 100)
#   for (i in 1:nrow(x1)) {
#     if (i %% 1000 == 0) print(i)
#     t.chrom <- x1[i, "chrom"] == as.character(x2[, "chrom"])
#     t.start <- (x1[i, "start"] > x2[, "start"])  & (x1[i, "start"] < x2[, "end"])
#     t.end <- (x1[i, "end"] > x2[, "start"])  & (x1[i, "end"] < x2[, "end"])
#     t.ids <- which(t.chrom & (t.start | t.end))
#     if (length(t.ids) != 0) y[i, 1:length(t.ids)] <- t.ids
#   }
#   return(y)
# }
