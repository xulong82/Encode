# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network
# Rev: August 18, 2014

# GENOTYPE, EPI-HOTSPOTS AND EPI-GENOTYPE DEFINITION

library(amap)
library(ape)
library(ggplot2)

#----------------------
rm(list = ls())
load("~/Dropbox/Network/R2/rna.rdt")  # rna
n.sp <- ncol(rna.f)  # sample size
ids.sp <- colnames(rna.f)  # samples
rna.z <- t(apply(rna.f, 1, scale))  # z-transform
colnames(rna.z) <- ids.sp
hc1 <- hcluster(t(rna.z), method = "pearson", link = "average")  # clust samples
hc2 <- hcluster(rna.z, method = "pearson", link = "average")  # clust genes
module <- cutree(hc2, h = .25)  # pho > .75
table(names(module) == rownames(rna.z)) 

pdf("~/Dropbox/Network/R2.fig/phylo1.pdf")
# plot(as.phylo(hc2))
plot(as.phylo(hc1), type = "unrooted", lab4ut = "axial")
dev.off()

ids.mod <- unique(module)
n.mod <- length(ids.mod)
rna.mod <- matrix(nrow = n.mod, ncol = n.sp, dimnames = list(ids.mod, ids.sp))
ids.rna.mod <- list()
for (i in 1:n.mod) {
  module1 <- module[module == i]
  ids.gene1 <- names(module1)
  data1 <- rna.z[ids.gene1, ]
  rna.mod[i, ] <- ifelse(rep(length(ids.gene1) == 1, n.sp), data1, colMeans(data1))
  ids.rna.mod[[i]] <- ids.gene1
}
save(rna.mod, file = "~/Dropbox/Network/R2/rna.mod.rdt")  # marks.bi
save(ids.rna.mod, file = "~/Dropbox/Network/R2/ids.rna.mod.rdt")  # marks.bi

#----------------------
load("~/Dropbox/Network/R2/mark1.rdt")  # marks.bi
source("~/Dropbox/Network/R2/analysis1x.R")

marks.bi <- marks.bi[, 1:4, ]  # take off TFs
n.mks <- dim(marks.bi)[2] # histone mark number: no TFs
n.genome <- dim(marks.bi)[3]  # genome length
marks <- dimnames(marks.bi)[[2]]

#--- define epigenomic hotspots
idx.mk1 <- matrix(nrow = n.genome, ncol = n.mks)  # new binary data per mark
for (i in 1:n.mks) {
  cat(i, "/", n.mks, marks[i], "\n")  # progress
  idx1 <- as.logical(colSums(marks.bi[, i, ]))
  idx.mk1[, i] <- as.numeric(idx1)  # new merged binary
}
idx <- as.logical(rowSums(idx.mk1))  # index
epi.mks <- myBed(which(idx))  # define the epi-geno regions
epi.mks <- epi.mks[(epi.mks$end - epi.mks$start) > 5, ]  # min: 1k
epi.mks <- epi.mks[(epi.mks$end - epi.mks$start) < 25, ]  # max: 10k
rownames(epi.mks) <- paste("epi", 1:nrow(epi.mks), sep = "")
save(epi.mks, file = "~/Dropbox/Network/R2/hotspots.rdt") 

colnames(idx.mk1) <- marks
hc1 <- hcluster(t(idx.mk1), method = "pearson", link = "average")  # clust samples
plot(as.phylo(hc1))
plot(as.phylo(hc1), type = "unrooted", lab4ut = "axial")

#--- FIGURE: density of the hotspots across the genome
epi1.dat <- data.frame(y = c(idx.mk1), mark = rep(marks, each = nrow(idx.mk1)))
epi1.dat$idx <- rep(1:nrow(idx.mk1), length(marks))
epi1.dat <- epi1.dat[epi1.dat$y == 1, ]
epi1.dat$pos <- epi1.dat$idx * 2e2 * 1e-6  # Mbp
pdf(file = "~/Dropbox/Network/R2.fig/epi1Hist.pdf", width = 12, height = 8)
ggplot(epi1.dat, aes(x = pos, colour = mark)) +
  geom_freqpoly(aes(group = mark), size = 2, binwidth = 20) +
# geom_density(size = 2, aes(color = mark, y = ..count..)) +
  theme_bw() +
  xlab("Chromosome") + ylab("Count") +
  scale_x_continuous(breaks = chrmid, labels = names(chrlen)) +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 18, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

#--- FIGURE: distribution of epigenomic-hotspots
epi2.dat <- data.frame(pos1 = 1:n.gnm, value = as.numeric(idx),
                       chrom = rep(0, n.gnm), pos2 = rep(0, n.gnm))

idx.chr <- rep(0, length(idx))
chrlen1 <- c(1, as.integer(chrlen * 1e6 %/% 200)) # 200bp unit
for (i in 1:length(chrlen)) idx.chr[chrlen1[i]:chrlen1[i+1]] <- i
epi.mks$chrom <- idx.chr[epi.mks$start]
epi.mks$pos <- (epi.mks$start + epi.mks$end) / 2 - chrlen1[epi.mks$chrom]
epi.mks$length <- epi.mks$end - epi.mks$start
epi.mks <- epi.mks[epi.mks$chrom != 21, ]

pdf(file = "~/Dropbox/Network/R2.fig/epiHotspot.pdf", width = 15, height = 10)
ggplot(epi.mks, aes(x = pos, y = 1)) +
  geom_point(aes(color = length, size = length)) +
# geom_point(color = "red", size = .1) +
# geom_vline(aes(xintercept = pos)) +
  facet_grid(chrom ~ .) +
  theme_bw() +
  xlab("Genome (200 bp unit)") + ylab("") +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.title.y = element_blank()) +
  theme(legend.position = "none", legend.direction = "horizontal", 
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"), legend.key = element_blank()) 
#       legend.title = element_blank(), legend.key = element_blank()) 
dev.off()
