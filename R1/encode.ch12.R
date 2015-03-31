# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - import the encode data
# Rev: May 14, 2014

library(parallel)
source("~/Dropbox/Encode/R/myfunction1.R")
basedir <- "~/Dropbox/Encode"

#--- ChIP-seq broadPeak data from Encode -----------------------------------------------------
setwd("~/Dropbox/Encode/ChIP-seq/broadPeak")

ch12.h3k4me1 <- mycolnames(read.delim("wgEncodePsuHistoneCh12H3k04me1FImmortal2a4bInputPk.broadPeak", header = F))
ch12.h3k4me3 <- mycolnames(read.delim("wgEncodePsuHistoneCh12H3k04me3FImmortal2a4bInputPk.broadPeak", header = F))
ch12.h3k9me3 <- mycolnames(read.delim("wgEncodePsuHistoneCh12H3k09me3FImmortal2a4bInputPk.broadPeak", header = F))
ch12.h3k27me3 <- mycolnames(read.delim("wgEncodePsuHistoneCh12H3k27me3FImmortal2a4bInputPk.broadPeak", header = F))
ch12.h3k36me3 <- mycolnames(read.delim("wgEncodePsuHistoneCh12H3k36me3FImmortal2a4bInputPk.broadPeak", header = F))

#--- Segment ENCODE broadPeak file -----------------------------------------------------------
setwd(basedir)

ch12.marks.bi <- cbind(myBinary(ch12.h3k4me1), myBinary(ch12.h3k4me3), myBinary(ch12.h3k9me3), 
                       myBinary(ch12.h3k27me3), myBinary(ch12.h3k36me3))
colnames(ch12.marks.bi) <- c("ch12-h3k4me1", "ch12-h3k4me3", "ch12-h3k9me3", "ch12-h3k27me3", "ch12-h3k36me3")
rownames(ch12.marks.bi) <- 1:nrow(ch12.marks.bi)

ch12.marks.si <- cbind(mySeg(ch12.h3k4me1), mySeg(ch12.h3k4me3), mySeg(ch12.h3k9me3), 
                       mySeg(ch12.h3k27me3), mySeg(ch12.h3k36me3))
colnames(ch12.marks.si) <- c("ch12-h3k4me1", "ch12-h3k4me3", "ch12-h3k9me3", "ch12-h3k27me3", "ch12-h3k36me3")
rownames(ch12.marks.si) <- 1:nrow(ch12.marks.si)

refgene.pm.bi <- myBinary(refgene.pm) 
refgene.tx.bi <- myBinary(refgene.tx)
refgene.ex.bi <- myBinary(refgene.ex)

#--- RNA-seq ----------------------------------------------------------------------------------
ens2sym <- read.delim("~/Dropbox/Genome/ensembl2symbol.map", header = F, row.names = 1)
ch12.1.rna <- read.delim("~/Dropbox/Encode/RSEM/CH12Rep1.genes.results", row.names = 1)[, -1]
ch12.2.rna <- read.delim("~/Dropbox/Encode/RSEM/CH12Rep2.genes.results", row.names = 1)[, -1]
ch12.rna <- (ch12.1.rna + ch12.2.rna) / 2  # Average
ch12.rna$symbol <- ens2sym[rownames(ch12.rna), ]
ch12.rna <- aggregate(. ~ symbol, ch12.rna, sum)[, -3]
ch12.rna <- data.frame(row.names = ch12.rna$symbol, ch12.rna[, 2:5])
rm(ch12.1.rna, ch12.2.rna, ens2sym)

# save(ch12.marks.bi, ch12.marks.si, refgene.pm.bi, refgene.tx.bi, refgene.ex.bi, ch12.rna, file = "~/Dropbox/Encode/R/ch12.RData")
# load("~/Dropbox/Encode/R/ch12.RData")

gene.bed <- read.delim("~/Dropbox/Info/symbol.bed")
gene.bed <- gene.bed[!duplicated(gene.bed$Associated.Gene.Name), ]
gene.bed <- mycolnames(data.frame(row.names = gene.bed$Associated.Gene.Name, gene.bed[, 2:4]))

gene.bed <- gene.bed[intersect(rownames(gene.bed), rownames(ch12.rna)), ]
gene.bed <- gene.bed[-grep("[GLJHM]", gene.bed$chrom), ]

ch12.rna <- ch12.rna[rownames(gene.bed), ]
