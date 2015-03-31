# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - HPC computing

library(parallel)
load("/data/xwang/R/tohpc.RData")

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
myParallel1 <- function(x) {  # Parallel run of myCount1
  # x: binary matrix
  # y: number of consecutive 1 in each x column
  time.start <- strptime(date(), "%a %b %d %H:%M:%S %Y")
  x.list <- list()
  for (i in 1:ncol(x)) x.list[[i]] <- x[, i]
  y <- mclapply(x.list, myCount1, mc.cores = 7)
  y <- matrix(unlist(y), nrow = 1)
  colnames(y) <- colnames(x)
  time.end <- strptime(date(), "%a %b %d %H:%M:%S %Y")
  difftime(time.end, time.start, units = "auto")
  return(y)
}
#----------------------------------------------------------------------------------------

marks.nb <- myParallel1(marks.bi)
marks.pm.nb <- myParallel1(marks.pm.bi)
marks.tx.nb <- myParallel1(marks.tx.bi)
marks.ex.nb <- myParallel1(marks.ex.bi)

save(marks.nb, marks.pm.nb, marks.tx.nb, marks.ex.nb, 
     file = "/data/xwang/R/outhpc.RData")

