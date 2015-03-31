# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - import the encode data
# Rev: May 2, 2014

library(parallel)
source("~/Dropbox/Encode/R/myfunction1.R")
basedir <- "~/Dropbox/Encode"

#--- ChIP-seq broadPeak data from Encode -----------------------------------------------------
setwd("~/Dropbox/Encode/ChIP-seq/broadPeak")

mel.h3k4me1 <- mycolnames(read.delim("wgEncodeLicrHistoneMelH3k04me1MImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k4me3 <- mycolnames(read.delim("wgEncodeLicrHistoneMelH3k04me3MImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k9ac <- mycolnames(read.delim("wgEncodeLicrHistoneMelH3k09acMImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k27ac <- mycolnames(read.delim("wgEncodeLicrHistoneMelH3k27acMImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k27me3 <- mycolnames(read.delim("wgEncodeLicrHistoneMelH3k27me3MImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k36me3 <- mycolnames(read.delim("wgEncodeLicrHistoneMelH3k36me3MImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k79me2 <- mycolnames(read.delim("wgEncodeLicrHistoneMelH3k79me2MImmortalC57bl6StdPk.broadPeak", header = F))

#--- Segment ENCODE broadPeak file -----------------------------------------------------------
setwd(basedir)

mel.marks.bi <- cbind(myBinary(mel.h3k4me1), myBinary(mel.h3k4me3), myBinary(mel.h3k9ac), 
                      myBinary(mel.h3k27ac), myBinary(mel.h3k27me3), myBinary(mel.h3k36me3), myBinary(mel.h3k79me2))
colnames(mel.marks.bi) <- c("MEL-h3k4me1", "MEL-h3k4me3", "MEL-h3k9ac", "MEL-h3k27ac", "MEL-h3k27me3", "MEL-h3k36me3", "MEL-h3k79me2")
rownames(mel.marks.bi) <- 1:nrow(mel.marks.bi)

mel.marks.si <- cbind(mySeg(mel.h3k4me1), mySeg(mel.h3k4me3), mySeg(mel.h3k9ac), 
                      mySeg(mel.h3k27ac), mySeg(mel.h3k27me3), mySeg(mel.h3k36me3), mySeg(mel.h3k79me2))
colnames(mel.marks.si) <- c("MEL-h3k4me1", "MEL-h3k4me3", "MEL-h3k9ac", "MEL-h3k27ac", "MEL-h3k27me3", "MEL-h3k36me3", "MEL-h3k79me2")
rownames(mel.marks.si) <- 1:nrow(mel.marks.si)

refgene.pm.bi <- myBinary(refgene.pm) 
refgene.tx.bi <- myBinary(refgene.tx)
refgene.ex.bi <- myBinary(refgene.ex)

#--- RNA-seq ----------------------------------------------------------------------------------
ens2sym <- read.delim("~/Dropbox/Genome/ensembl2symbol.map", header = F, row.names = 1)
# htseq <- read.delim("~/Dropbox/Encode/HTSEQ/wgEncodeLicrRnaSeqMelCellPapMImmortalC57bl6AlnRep1.txt", header = F)
mel1.rna <- read.delim("~/Dropbox/Encode/RSEM/MELRep1_cut_0.genes.results", row.names = 1)[, -1]
mel2.rna <- read.delim("~/Dropbox/Encode/RSEM/MELRep2_cut_0.genes.results", row.names = 1)[, -1]
mel.rna <- (mel1.rna + mel2.rna) / 2  # Average
mel.rna$symbol <- ens2sym[rownames(mel.rna), ]
mel.rna <- aggregate(. ~ symbol, mel.rna, sum)[, -3]
mel.rna <- data.frame(row.names = mel.rna$symbol, mel.rna[, 2:5])
rm(mel1.rna, mel2.rna, ens2sym)

# save(mel.marks.bi, mel.marks.si, refgene.pm.bi, refgene.tx.bi, refgene.ex.bi, mel.rna, file = "~/Dropbox/Encode/R/mel.RData")
# load("~/Dropbox/Encode/R/mel.RData")

gene.bed <- read.delim("~/Dropbox/Genome/symbol.bed")
gene.bed <- gene.bed[!duplicated(gene.bed$Associated.Gene.Name), ]
gene.bed <- mycolnames(data.frame(row.names = gene.bed$Associated.Gene.Name, gene.bed[, 2:4]))

gene.bed <- gene.bed[intersect(rownames(gene.bed), rownames(mel.rna)), ]
gene.bed <- gene.bed[-grep("[GLJHM]", gene.bed$chrom), ]

mel.rna <- mel.rna[rownames(gene.bed), ]
