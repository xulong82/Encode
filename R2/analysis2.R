# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenetic network
# Rev: August 18, 2014

library(ggplot2)

rm(list = ls())
load("~/Dropbox/Network/R2/scan1f.rdt")
load("~/Dropbox/Network/R2/rna.mod.rdt")  # rna
load("~/Dropbox/Network/R2/ids.rna.mod.rdt")  # mod:gene
load("~/Dropbox/Network/R2/hotspots.rdt")  # epi-hotspots
load("~/Dropbox/Network/R2/epigeno.rdt")  # epi-genotype
source("~/Dropbox/Network/R2/analysis2x.R")
#---------------------------------------------------------------------------
n.mod <- nrow(rna.mod)
id.mod <- rownames(rna.mod)

scan1.h3k4me1 <- matrix(nrow = 0, ncol = 6, 
                        dimnames = list(NULL, c("mod", "pho", "p", "q", "beta", "logLik")))
scan1.h3k27me3 <- scan1.h3k27ac <- scan1.h3k4me3 <- scan1.h3k4me1
scan1.h3k4me1.2 <- matrix(nrow = 0, ncol = 7, 
                        dimnames = list(NULL, c("gene", "mod", "pho", "p", "q", "beta", "logLik")))
scan1.h3k27me3.2 <- scan1.h3k27ac.2 <- scan1.h3k4me3.2 <- scan1.h3k4me1.2
scan1.mmark <- matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c("mod", "epi")))
scan1.mmark.2 <- matrix(nrow = 0, ncol = 2, dimnames = list(NULL, c("gene", "epi")))
for (i in 1:length(scan1f)) {
  if (i %% 1e1 == 0) cat(i, "in", length(scan1f), "\n")
  gene1 <- ids.rna.mod[[i]]
  
  if (length(dim(scan1f[[i]])) == 3) {
    mmark1 <- as.data.frame(t(scan1f[[i]][, , "q"]))
    if (nrow(mmark1) > 0) {
      for (j in 1:nrow(mmark1)) {
        mmark11 <- sort(mmark1[j, ])
        idx.epi <- rownames(mmark1)
        if (length(mmark11) > 1) {
          if (mmark11[2] < .05) {
            mmark11.1 <- data.frame(mod = i, epi = idx.epi[j])
            mmark11.2 <- data.frame(gene = gene1, epi = rep(idx.epi[j], length(gene1)))
            scan1.mmark <- rbind(scan1.mmark, mmark11.1)
            scan1.mmark.2 <- rbind(scan1.mmark.2, mmark11.2)
          }
        }
      }
    }
    # --- code for 1 mark 
  }
}
save(scan1.mmark, scan1.mmark.2, file = "~/Dropbox/Network/R2/mmark.rdt")
#     mark1 <- as.data.frame(scan1f[[i]]["h3k4me1", , ])
#     mark1 <- mark1[!is.na(mark1[, "q"]), ]
#     mark1 <- mark1[mark1[, "q"] < .05, ]
#     mark1 <- cbind(mod = rep(i, nrow(mark1)), mark1)
#     scan1.h3k4me1 <- rbind(scan1.h3k4me1, mark1)
#     mark2 <- data.frame(gene = rep(gene1, each = nrow(mark1)), 
#                         epi = rep(rownames(mark1), length(gene1)), 
#                         logLik = rep(mark1$logLik, length(gene1)))
#     scan1.h3k4me1.2 <- rbind(scan1.h3k4me1.2, mark2)
#     
#     mark1 <- as.data.frame(scan1f[[i]]["h3k4me3", , ])
#     mark1 <- mark1[!is.na(mark1[, "q"]), ]
#     mark1 <- mark1[mark1[, "q"] < .05, ]
#     mark1 <- cbind(mod = rep(i, nrow(mark1)), mark1)
#     scan1.h3k4me3 <- rbind(scan1.h3k4me3, mark1)
#     mark2 <- data.frame(gene = rep(gene1, each = nrow(mark1)), 
#                         epi = rep(rownames(mark1), length(gene1)), 
#                         logLik = rep(mark1$logLik, length(gene1)))
#     scan1.h3k4me3.2 <- rbind(scan1.h3k4me3.2, mark2)
#     
#     mark1 <- as.data.frame(scan1f[[i]]["h3k27ac", , ])
#     mark1 <- mark1[!is.na(mark1[, "q"]), ]
#     mark1 <- mark1[mark1[, "q"] < .05, ]
#     mark1 <- cbind(mod = rep(i, nrow(mark1)), mark1)
#     scan1.h3k27ac <- rbind(scan1.h3k27ac, mark1)
#     mark2 <- data.frame(gene = rep(gene1, each = nrow(mark1)), 
#                         epi = rep(rownames(mark1), length(gene1)), 
#                         logLik = rep(mark1$logLik, length(gene1)))
#     scan1.h3k27ac.2 <- rbind(scan1.h3k27ac.2, mark2)
#     
#     mark1 <- as.data.frame(scan1f[[i]]["h3k27me3", , ])
#     mark1 <- mark1[!is.na(mark1[, "q"]), ]
#     mark1 <- mark1[mark1[, "q"] < .05, ]
#     mark1 <- cbind(mod = rep(i, nrow(mark1)), mark1)
#     scan1.h3k27me3 <- rbind(scan1.h3k27me3, mark1)
#     mark2 <- data.frame(gene = rep(gene1, each = nrow(mark1)), 
#                         epi = rep(rownames(mark1), length(gene1)), 
#                         logLik = rep(mark1$logLik, length(gene1)))
#     scan1.h3k27me3.2 <- rbind(scan1.h3k27me3.2, mark2)

#---------------------------------------------------------------------------
mark <- "any"
mark <- "mmark"

id.epi1 <- idx.epi[[mark]][[1]]
len.epi <- length(id.epi1)
id.gene1 <- ids.rna.mod[[1]]
dat1 <- data.frame(mod = rep(id.mod[1], length(id.epi1)), epi = id.epi1)
dat1$mod <- as.character(dat1$mod)
map1.dat <- data.frame(gene = rep(id.gene1, each = nrow(dat1)), 
                       mod = rep(dat1$mod, length(id.gene1)), 
                       epi = rep(dat1$epi, length(id.gene1))) 
for (i in 2:n.mod) {  # per mod
  if (i %% 50 == 0) cat(i, "in", n.mod, "\n")
  id.epi1 <- idx.epi[[mark]][[i]]
  len.epi <- c(len.epi, length(id.epi1))
  id.gene1 <- ids.rna.mod[[i]]
  dat1 <- data.frame(mod = rep(id.mod[i], length(id.epi1)), epi = id.epi1)
  dat1$mod <- as.character(dat1$mod)
  dat2 <- data.frame(gene = rep(id.gene1, each = nrow(dat1)), 
                     mod = rep(dat1$mod, length(id.gene1)), 
                     epi = rep(dat1$epi, length(id.gene1))) 
  map1.dat <- rbind(map1.dat, dat2)
}

map1.dat$gene <- as.character(map1.dat$gene)
map1.dat$epi.pos <- rowMeans(epi.mks[map1.dat$epi, ]) * 2e2 * 1e-6  # in Mbp
map1.dat$gene.pos <- ens.ucsc[map1.dat$gene, "pos.Mbp"]

hist(len.epi, n = length(len.epi))
table(len.epi > 0)
len.epi[len.epi > 0]
idx.epi[[mark]][len.epi > 10]

which(len.epi == max(len.epi))
idx1 <- idx.epi[[mark]][[which(len.epi == max(len.epi))]]
idx1.epi <- epi.mks[idx1, ]
summary(idx1.epi$end - idx1.epi$start)
hist(idx1.epi$end - idx1.epi$start)

idx1.epigeno <- epigeno[, , idx1]
idx1.epigeno[, 1, 1:15]
cor(idx1.epigeno[, 1, 1:15])

#---------------------------------------------------------------------------
mark <- "h3k4me1"
mark <- "h3k4me3"
mark <- "h3k27ac"
mark <- "h3k27me3"

id.epi1 <- idx.epi$pmark[[1]][[mark]]
len.epi <- length(id.epi1)
id.gene1 <- ids.rna.mod[[1]]
dat1 <- data.frame(mod = rep(id.mod[1], length(id.epi1)), epi = id.epi1)
dat1$mod <- as.character(dat1$mod)
map1.dat <- data.frame(gene = rep(id.gene1, each = nrow(dat1)), 
                       mod = rep(dat1$mod, length(id.gene1)), 
                       epi = rep(dat1$epi, length(id.gene1))) 
for (i in 2:n.mod) {  # per mod
  if (i %% 50 == 0) cat(i, "in", n.mod, "\n")
  if (length(idx.epi$pmark[[i]]) != 0) {
    id.epi1 <- idx.epi$pmark[[i]][[mark]]
    len.epi <- c(len.epi, length(id.epi1))
    id.gene1 <- ids.rna.mod[[i]]
    dat1 <- data.frame(mod = rep(id.mod[i], length(id.epi1)), epi = id.epi1)
    dat1$mod <- as.character(dat1$mod)
    dat2 <- data.frame(gene = rep(id.gene1, each = nrow(dat1)), 
                       mod = rep(dat1$mod, length(id.gene1)), 
                       epi = rep(dat1$epi, length(id.gene1))) 
    map1.dat <- rbind(map1.dat, dat2)
  } else {
    len.epi <- c(len.epi, 0)
  }
}
map1.dat$gene <- as.character(map1.dat$gene)
map1.dat$epi.pos <- rowMeans(epi.mks[map1.dat$epi, ]) * 2e2 * 1e-6  # in Mbp
map1.dat$gene.pos <- ens.ucsc[map1.dat$gene, "pos.Mbp"]

hist(len.epi, n = length(len.epi))
table(len.epi > 0)
which(len.epi == max(len.epi))
idx1 <- idx.epi$pmark[[which(len.epi == max(len.epi))]][[mark]]

idx1.epi <- epi.mks[idx1, ]
summary(idx1.epi$end - idx1.epi$start)
hist(idx1.epi$end - idx1.epi$start)

idx1.epigeno <- epigeno[, , idx1]

cor(idx1.epigeno[, 4, 1:10])

cutoff.pho <- .9
cutoff.pval  <- .01
cutoff.beta <- .5
cutoff.lod <- 3
idx.epi <- list()

#--- per gene: epi-spots with required p and pho on 1 or more marks
# for (i in 1:length(cortest)) {
#   if (i %% 1e2 == 0) cat(i, "in", length(cortest), "\n")
#   temp1 <- cortest[[i]]
#   temp.p <- temp1[, , "p"]  # p value
#   temp.pho <- temp1[, , "pho"]  # pho value
#   null <- which(apply(temp.p, 2, function (x) {length(which(is.na(x))) == 4}))
#   if (length(null) != 0) {
#     temp.p <- temp.p[, -null]  # zero-variation in epigenotype
#     temp.pho <- temp.pho[, -null]
#   }
#   #--- significant spots in any mark
#   true.p <- apply(temp.p, 2, function (x) {min(x, na.rm = T) < cutoff.pval})
#   true.pho <- apply(temp.pho, 2, function (x) {max(abs(x), na.rm = T) > cutoff.pho})
#   idx.epi$any[[i]] <- unname(which(true.p & true.pho))
#   #--- significant spots per mark
#   true.p <- apply(temp.p, 1, function(x) {x < cutoff.pval})
#   true.pho <- apply(temp.pho, 1, function(x) {abs(x) > cutoff.pho})
#   idx.epi$pmark[[i]] <- apply(true.p & true.pho, 2, function (x) {unname(which(x))})
#   #--- spots with 2 or more significant marks
#   idx.mmark <- 0
#   for (j in 1:ncol(temp.p)) {
#     temp.p1 <- sort(temp.p[, j])
#     temp.pho1 <- sort(abs(temp.pho[, j]), decreasing = T)
#     if (length(temp.p1) > 1) {
#       if ((temp.p1[2] < cutoff.pval) & (temp.pho1[2] > cutoff.pho)) 
#         idx.mmark <- c(idx.mmark, j)
#     }
#   }
#   idx.epi$mmark[[i]] <- idx.mmark[-1]
# }
# cat(date(), "fi \n")
# save(idx.epi, file = "/data/xwang/Network/R/scan1/idx.epi.rdt")
# 
