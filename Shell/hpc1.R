# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network
# Rev: August 18, 2014

cat(date(), "DEFINE EPIGENOTYPE VALUES PER MARK PER SAMPLE \n")

library(parallel)

load("/data/xwang/Network/R/mark1.rdt")
load("/data/xwang/Network/R/hotspots.rdt")
#------------------------------------
marks.bi <- marks.bi[, 1:4, ]  # histone marks only
n.sps <- dim(marks.bi)[1]  # sample number
n.mks <- dim(marks.bi)[2]  # mark number
n.epi <- nrow(epi.mks)  # epispot number
id.sps <- dimnames(marks.bi)[[1]]  # sample name
id.mks <- dimnames(marks.bi)[[2]]  # mark name
id.epi <- rownames(epi.mks)  # epispot name
epi.ln <- epi.mks$end - epi.mks$start  # hotspot length
#------------------------------------
myEpigeno <- function (x) {  # per sample
  epigeno <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, id.mks))
  for (i in 1:n.mks) {
    cat(i, "/", n.mks, id.mks[i], "\n")  # marks
    marks.bi1 <- x[i, ]
    for (j in 1:n.epi) {
      if (j %% 1e4 == 0) cat(j, "/", n.epi, "\n")
      epigeno[j, i] <- sum(marks.bi1[epi.mks[j, 1]:epi.mks[j, 2]])
    }
  }
  epigeno <- epigeno / epi.ln
  return(epigeno)
}
#------------------------------------
marks.per.sample <- list()
for (i in 1:n.sps) marks.per.sample[[i]] = marks.bi[i, , ]
y <- mclapply(marks.per.sample, myEpigeno, mc.cores = 11)
epigeno <- array(NA, c(n.sps, n.mks, n.epi), dimnames = list(id.sps, id.mks, id.epi))
for (i in 1:length(y)) epigeno[i, , ] = t(y[[i]])
#------------------------------------
marked <- rep(0, n.epi)
for (i in 1:n.epi) marked[i] <- length(which(as.logical(epigeno[, , i])))
epigeno <- epigeno[, , marked > 4]  # > 4 in 44 spots
save(epigeno, file = "/data/xwang/Network/R/epigeno.rdt")
cat(date(), "fi \n")

