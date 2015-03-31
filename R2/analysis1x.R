# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network - Regression 1
# Rev: June 7, 2014

#-------------------------------------------------
chr.ucsc <- read.delim("~/Dropbox/X/chromInfo.txt", header = F)  # chromosomes length
chrlen <- cumsum(as.numeric(chr.ucsc$V2)) * 1e-6
names(chrlen) <- 1:21
chrmid <- diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])

#-------
myBed <- function(x) {  # define the epigenome hotspots
  # x: index; y: bed data declares epigenome hotspots
  leng <- length(x)
  id <- 1  # id of region
  start <- x[1]  # the start index
  y0 <- matrix(nrow = 1e6, ncol = 2)  # start, end
  for (i in 2:leng) {
    if (i %% 1e5 == 0) cat(paste(as.integer(i * 100 / leng), "%", sep = ""), "\n")
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
  y0 <- y0[!is.na(y0[, 1]), ]  # truncation
  y1 <- data.frame(row.names = 1:nrow(y0), start = y0[, 1], end = y0[, 2])
  return(y1)
}

