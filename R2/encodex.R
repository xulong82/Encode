# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network - import the encode data
# Rev: May 30, 2014

#--- UCSC features
chrom.info <- read.delim("~/Dropbox/X/chromInfo.txt", header = F)  # chromosome length
chrom.unit <- 200  # unit: nucleotide length
chrom.leng <- chrom.info$V2 %/% chrom.unit + 1  # chromosome length in chrom.unit
chrom.bina <- rep(0, sum(chrom.leng))

refgene <- read.delim("~/Dropbox/X/myRefGene.tx", header = T, stringsAsFactors = F)
refgene <- refgene[, c("name2", "chrom", "txStart", "txEnd")]  # RefSeq transcripts
refgene <- refgene[!duplicated(refgene$name2), ]
refexon <- read.delim("myRefGene.exons", header = F)  # RefSeq exons
refgene.pm <- data.frame(row.names = refgene$name2,
  chrom = refgene$chrom,  # Promoter: -2k -> +2k TSS 
  chromStart = refgene$txStart - 2000,
  chromEnd   = refgene$txStart + 2000)
refgene.tx <- data.frame(row.names = refgene$name2,
  chrom = refgene$chrom,  # Promoter: -2k -> +2k TSS 
  chromStart = refgene$txStart,
  chromEnd   = refgene$txEnd)

#----------------------------------------------------------------------------------------
myColname <- function (x) {  # Assign colnames for ENCODE broadPeak file
  colnames(x)[1:3] <- c("chrom", "chromStart", "chromEnd") 
  return(x)
}
#----------------------------------------------------------------------------------------
myBinary <- function(x) {  # Generate binary file for ENCODE broadPeak file
  # x: BED columns 1-3 (1:chrom; 2:chromStart; 3:chromEnd)
  # y: binary signal of BED ranges
  y <- chrom.bina  # Initialize the return genome
  x$chrom <- gsub("chrX", "chr20", x$chrom)
  x$chrom <- gsub("chrY", "chr21", x$chrom)  # Encode has no Y chromosome
  x$chrom <- as.numeric(gsub("chr", "", x$chrom))
  for (i in 1:nrow(x)) {
    if (i %% 1e3 == 0) cat(i, "in", nrow(x), "\n")
    x1 <- x[i, ]  # Single record
    offset <- ifelse(x1$chrom == 1, 0, sum(chrom.leng[1:(x1$chrom - 1)]))
    start <- offset + x1$chromStart %/% chrom.unit
    end <- offset + x1$chromEnd %/% chrom.unit
    y[start:end] <- 1
  }
  return(y)  # genome wide binary signal of the input feature
}
#----------------------------------------------------------------------------------------
mySignal <- function(x) {  # Generate signal intensity file from ENCODE broadPeak file
  # x: ENCODE broadPeak file
  # y: segmentized broadPeak file with 200bp unit
  y <- chrom.bina  # Initialize the return genome
  x$chrom <- gsub("chrX", "chr20", x$chrom)
  x$chrom <- gsub("chrY", "chr21", x$chrom)  # Encode has no Y chromosome
  x$chrom <- as.numeric(gsub("chr", "", x$chrom))
  for (i in 1:nrow(x)) {
    if (i %% 1e3 == 0) cat(i, "in", nrow(x), "\n")
    x1 <- x[i, ]  # Single record
    offset <- ifelse(x1$chrom == 1, 0, sum(chrom.leng[1:(x1$chrom - 1)]))
    start <- offset + x1$chromStart %/% chrom.unit
    end <- offset + x1$chromEnd %/% chrom.unit
    y[start:end] <- x1$V7
  }
  return(y)  # genome wide binary signal of the input feature
}
#----------------------------------------------------------------------------------------
