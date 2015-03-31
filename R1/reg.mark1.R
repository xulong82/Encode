# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenetic network
# Rev: May 14, 2014
#---------------------------------------------------------------------------
library(MASS)
library(ggplot2)
#---------------------------------------------------------------------------
load("~/Dropbox/Encode/R/marks.rdt")
load("~/Dropbox/Encode/R/marks.broadpeak.rdt")
load("~/Dropbox/Encode/R/refgene.rdt")
source("~/Dropbox/Encode/R/myfunction1.R")
#----------------------
epi.mks = list()  # define epi-geno per mark: region
for (i in 1:n.mks) {
  cat(i, "/", n.mks, marks[i], "\n")  # progress
  idx.mks <- which(as.logical(colSums(marks.bi[, i, ])))
  epi.mks[[i]] <- myEpigeno(idx.mks)
  rownames(epi.mks[[i]]) <- paste(marks[i], rownames(epi.mks[[i]]), sep = ":")
  epi.mks[[i]] <- epi.mks[[i]][epi.mks[[i]]$start != epi.mks[[i]]$end, ]  # delete length-1 entries
}
names(epi.mks) <- marks
#----------------------
epi.sps <- list()  # assign epi-genotype for samples 
for (i in 1:n.mks) {
  epi.mks0 <- epi.mks[[i]]
  marks.bi0 <- marks.bi[, i, ]
# marks.si0 <- marks.si[, i, ]
  n.epi <- nrow(epi.mks0)
  epi.sps0 <- matrix(nrow = n.epi, ncol = n.sps, dimnames = list(rownames(epi.mks0), samples))
  for (j in 1:n.sps) {
    cat(j, "/", n.sps, samples[j], "\n")  # progress
    for (k in 1:n.epi) {
      if (k %% 1000 == 0) cat(k, "/", n.epi, "\n")
      start <- epi.mks0[k, 1]
      end <- epi.mks0[k, 2]
      epi.sps0[k, j] <- sum(marks.bi0[j, start:end])
#     epi.sps0[k, j] <- mean(marks.si0[j, start:end])
    }
  }
  epi.sps[[i]] <- epi.sps0
}
names(epi.sps) <- marks
save(epi.mks, epi.sps, file = "~/Dropbox/Encode/R/epi1.rdt")
#---------------------------------------------------------------------------
