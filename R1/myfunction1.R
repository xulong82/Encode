# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - define the functions
# Rev: March 13, 2014

source("~/Dropbox/Encode/R/myinfo.R")
#----------------------------------------------------------------------------------------
mycolnames <- function (x) {  # Assign colnames for ENCODE broadPeak file
  colnames(x)[1:3] <- c("chrom", "chromStart", "chromEnd") 
  return(x)
}
#----------------------------------------------------------------------------------------
myBinary <- function(x) {  # Generate binary file from ENCODE broadPeak file
  # x: BED columns 1-3 (1:chrom; 2:chromStart; 3:chromEnd)
  # y: binary signal of BED ranges
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
    y[start:end] <- 1
  }
  close(pb)  # Close progress bar
  return(y)  # genome wide binary signal of the input feature
}
#----------------------------------------------------------------------------------------
mySignal <- function(x) {  # Generate signal intensity file from ENCODE broadPeak file
  # x: ENCODE broadPeak file
  # y: segmentized broadPeak file with 200bp unit
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
myEpigeno <- function(x) {  # Define the epi-genotypes bed data 
  # x: binary sequence; y: bed data declares epi-geno regions
  id <- 1  # epi-geno ID
  start <- x[1]  # start index
  y0 <- matrix(nrow = 1000000, ncol = 2)  # start, end
  pb <- txtProgressBar(min = 0, max = 100, style = 3)  # Progress bar
  for (i in 2:length(x)) {
    if (i %% 10000 == 0) {
      progress <- (i * 100) %/% length(x)
      setTxtProgressBar(pb, progress)  
    }  # Update progress bar
    diff <- x[i] - x[i-1]
    if (diff != 1) {
      end <- x[i-1]
      y0[id, ] <- c(start, end)
      start <- x[i]
      id <- id + 1
    } 
  }
  end <- x[i]
  y0[id, ] <- c(start, end)
  close(pb)  # Close progress bar
  y0 <- y0[!is.na(y0[, 1]), ]  # truncation
  y1 <- data.frame(row.names = 1:nrow(y0), start = y0[, 1], end = y0[, 2])
  return(y1)
}
#----------------------------------------------------------------------------------------
myLength <- function(x) {  # Calculate length of consecutive 1 sequences
  # x: binary sequence; y: length of consecutive 1 sequences
  y0 <- ifelse(x[1] == 0, 0, 1)  # status: 1
  id <- ifelse(x[1] == 0, 0, 1)  # id for sequences
  y1 <- rep(NA, 100000)  # initiate vecotor for length of 1 sequences
  pb <- txtProgressBar(min = 0, max = 100, style = 3)  # Progress bar
  len <- length(x)
  for (i in 2:len) {
    if (i %% 100000 == 0) {
      progress <- (i * 100) %/% length(x)
      setTxtProgressBar(pb, progress)  
    }  # Update progress bar
    turn <- x[i] - x[i-1]
    if (y0) {
      if (turn == 0) {
        y0 = y0 + 1
        if (i == len) y1[id] = y0
      } else {
        y1[id] = y0
        y0 = 0
      } 
    } else {
      if (turn == 1) {
        y0 = 1
        id = id + 1
        if (i == len) y1[id] = y0
      }
      y0 <- ifelse(turn == 0, 0, 1)
    }
  }
  close(pb)  # Close progress bar
  return(y1)
}
#----------------------------------------------------------------------------------------
myCount1 <- function(x) {  # Calculate number of consecutive 1 sequences
  # x: binary sequence; y: number of consecutive 1 sequences
  y <- ifelse(x[1] == 0, 0, 1)  # initiate y
  pb <- txtProgressBar(min = 0, max = 100, style = 3)  # Progress bar
  for (i in 2:length(x)) {
    if (i %% 100000 == 0) {
      progress <- (i * 100) %/% length(x)
      setTxtProgressBar(pb, progress)  
    }  # Update progress bar
    turn <- x[i] - x[i-1]
    y <- y + ifelse(turn == 1, 1, 0)
  }
  close(pb)  # Close progress bar
  return(y)
}
#----------------------------------------------------------------------------------------
myCount2 <- function(x1, x2) {  # Calculate number of 1 or sum of signal in a region
  # x1: BED columns 1-3 (1:chrom; 2:chromStart; 3:chromEnd)
  # x2: Binary genome/epigenome features
  # y:  number of 1 or sum of signal in x2 with regions defined by x1
  y <- c(rep(0, length(x1)))
  x1$chrom <- gsub("chrX", "chr20", x1$chrom)
  x1$chrom <- gsub("chrY", "chr21", x1$chrom)  # Encode has no Y chromosome
  x1$chrom <- gsub("X", "20", x1$chrom)
  x1$chrom <- gsub("Y", "21", x1$chrom)  # Encode has no Y chromosome
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
myPattern <- function(x1, x2) {  # Assign pattern of marks to genes
  # x1: epigenome features
  # x2: pattern matrix
  y <- matrix(0, nrow = nrow(x1), ncol = nrow(x2))  # pattern of marks to genes
  for (i in 1:nrow(x1)) {
    if (i %% 1000 == 0) print(i)
    x1.1 <- x1[i, ]  # Single record
    x1.1[x1.1 > 1] <- 1
    index <- which(apply(x2, 1, function (x) {paste(x, collapse = "") == paste(x1.1, collapse = "")}))
    y[i, index] <- 1  # flexible
  }
  return(y)
}
#----------------------------------------------------------------------------------------
# refgene.pm.bi <- myBinary(refgene.pm) 
# refgene.tx.bi <- myBinary(refgene.tx)
# refgene.ex.bi <- myBinary(refgene.ex)
# save(refgene.pm.bi, refgene.tx.bi, refgene.ex.bi, file = "~/Dropbox/Encode/R/refgene.rdt")
