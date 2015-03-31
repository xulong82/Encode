# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network - import the encode data
# Rev: May 30, 2014

source("~/Dropbox/Network/R2/encodex.R")

#--- Encode: ChIP-seq broadPeak files
setwd("~/Dropbox/Network/ChIP-seq/marks")
cbellum.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneCbellumH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F))
cbellum.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneCbellumH3k27me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
cbellum.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneCbellumH3k4me1MAdult8wksC57bl6StdPk.broadPeak", header = F))
cbellum.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneCbellumH3k4me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
esb4.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneEsb4H3k27acME0C57bl6StdPk.broadPeak", header = F))
esb4.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneEsb4H3k27me3ME0C57bl6StdPk.broadPeak", header = F))
esb4.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneEsb4H3k4me1ME0C57bl6StdPk.broadPeak", header = F))
esb4.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneEsb4H3k4me3ME0C57bl6StdPk.broadPeak", header = F))
heart.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneHeartH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F))
heart.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneHeartH3k27me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
heart.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneHeartH3k4me1MAdult8wksC57bl6StdPk.broadPeak", header = F))
heart.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneHeartH3k4me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
kidney.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneKidneyH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F))
kidney.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneKidneyH3k27me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
kidney.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneKidneyH3k4me1MAdult8wksC57bl6StdPk.broadPeak", header = F))
kidney.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneKidneyH3k4me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
liver.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneLiverH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F))
liver.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneLiverH3k27me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
liver.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneLiverH3k4me1MAdult8wksC57bl6StdPk.broadPeak", header = F))
liver.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneLiverH3k4me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
mel.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneMelH3k04me1MImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneMelH3k04me3MImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneMelH3k27acMImmortalC57bl6StdPk.broadPeak", header = F))
mel.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneMelH3k27me3MImmortalC57bl6StdPk.broadPeak", header = F))
smint.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneSmintH3k04me1MAdult8wksC57bl6StdPk.broadPeak", header = F))
smint.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneSmintH3k04me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
smint.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneSmintH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F))
smint.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneSmintH3k27me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
spleen.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneSpleenH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F))
spleen.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneSpleenH3k27me3MAdlt8wC57bl6StdPk.broadPeak", header = F))
spleen.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneSpleenH3k4me1MAdult8wksC57bl6StdPk.broadPeak", header = F))
spleen.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneSpleenH3k4me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
testis.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneTestisH3k04me1MAdult8wksC57bl6StdPk.broadPeak", header = F))
testis.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneTestisH3k04me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
testis.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneTestisH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F))
testis.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneTestisH3k27me3MAdlt8wC57bl6StdPk.broadPeak", header = F))
thymus.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneThymusH3k04me1MAdult8wksC57bl6StdPk.broadPeak", header = F))
thymus.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneThymusH3k04me3MAdult8wksC57bl6StdPk.broadPeak", header = F))
thymus.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneThymusH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F))
thymus.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneThymusH3k27me3MAdlt8wC57bl6StdPk.broadPeak", header = F))
brain.h3k4me1 <- myColname(read.delim("wgEncodeLicrHistoneWbrainH3k04me1UE14halfC57bl6StdPk.broadPeak", header = F))
brain.h3k4me3 <- myColname(read.delim("wgEncodeLicrHistoneWbrainH3k04me3UE14halfC57bl6StdPk.broadPeak", header = F))
brain.h3k27ac <- myColname(read.delim("wgEncodeLicrHistoneWbrainH3k27acUE14halfC57bl6StdPk.broadPeak", header = F))
brain.h3k27me3 <- myColname(read.delim("wgEncodeLicrHistoneWbrainH3k27me3ME14halfC57bl6StdPk.broadPeak", header = F))
cbellum.ctcf <- myColname(read.delim("wgEncodeLicrTfbsCbellumCtcfMAdult8wksC57bl6StdPk.broadPeak", header = F))
cbellum.pol2 <- myColname(read.delim("wgEncodeLicrTfbsCbellumPol2MAdult8wksC57bl6StdPk.broadPeak", header = F))
esb4.ctcf <- myColname(read.delim("wgEncodeLicrTfbsEsb4CtcfME0C57bl6StdPk.broadPeak", header = F))
esb4.pol2 <- myColname(read.delim("wgEncodeLicrTfbsEsb4Pol2ME0C57bl6StdPk.broadPeak", header = F))
heart.ctcf <- myColname(read.delim("wgEncodeLicrTfbsHeartCtcfMAdult8wksC57bl6StdPk.broadPeak", header = F))
heart.pol2 <- myColname(read.delim("wgEncodeLicrTfbsHeartPol2MAdult8wksC57bl6StdPk.broadPeak", header = F))
kidney.ctcf <- myColname(read.delim("wgEncodeLicrTfbsKidneyCtcfMAdult8wksC57bl6StdPk.broadPeak", header = F))
kidney.pol2 <- myColname(read.delim("wgEncodeLicrTfbsKidneyPol2MAdult8wksC57bl6StdPk.broadPeak", header = F))
liver.ctcf <- myColname(read.delim("wgEncodeLicrTfbsLiverCtcfMAdult8wksC57bl6StdPk.broadPeak", header = F))
liver.pol2 <- myColname(read.delim("wgEncodeLicrTfbsLiverPol2MAdult8wksC57bl6StdPk.broadPeak", header = F))
mel.ctcf <- myColname(read.delim("wgEncodeLicrTfbsMelCtcfMImmortalC57bl6StdPk.broadPeak", header = F))
mel.pol2 <- myColname(read.delim("wgEncodeLicrTfbsMelPol2MImmortalC57bl6StdPk.broadPeak", header = F))
smint.ctcf <- myColname(read.delim("wgEncodeLicrTfbsSmintCtcfMAdult8wksC57bl6StdPk.broadPeak", header = F))
smint.pol2 <- myColname(read.delim("wgEncodeLicrTfbsSmintPol2MAdult8wksC57bl6StdPk.broadPeak", header = F))
spleen.ctcf <- myColname(read.delim("wgEncodeLicrTfbsSpleenCtcfMAdult8wksC57bl6StdPk.broadPeak", header = F))
spleen.pol2 <- myColname(read.delim("wgEncodeLicrTfbsSpleenPol2MAdult8wksC57bl6StdPk.broadPeak", header = F))
testis.ctcf <- myColname(read.delim("wgEncodeLicrTfbsTestisCtcfMAdult8wksC57bl6StdPk.broadPeak", header = F))
testis.pol2 <- myColname(read.delim("wgEncodeLicrTfbsTestisPol2MAdult8wksC57bl6StdPk.broadPeak", header = F))
thymus.ctcf <- myColname(read.delim("wgEncodeLicrTfbsThymusCtcfMAdult8wksC57bl6StdPk.broadPeak", header = F))
thymus.pol2 <- myColname(read.delim("wgEncodeLicrTfbsThymusPol2MAdult8wksC57bl6StdPk.broadPeak", header = F))
brain.ctcf <- myColname(read.delim("wgEncodeLicrTfbsWbrainCtcfUE14halfC57bl6StdPk.broadPeak", header = F))
brain.pol2 <- myColname(read.delim("wgEncodeLicrTfbsWbrainPol2UE14halfC57bl6StdPk.broadPeak", header = F))
save.image("~/Dropbox/Encode/R2/en.broadpeak.rdt")
load("~/Dropbox/Network/R2/broadpeak.rdt")

#--- ENCODE broadPeak file segments
samples <- c("cbellum", "esb4", "heart", "kidney", "liver", "mel", "smint", "spleen", "testis", "thymus", "brain")
marks <- c("h3k4me1", "h3k4me3", "h3k27ac", "h3k27me3", "ctcf", "pol2")
n.sps <- length(samples)
n.mks <- length(marks)
g.len <- sum(chrom.leng)

#--- I. Binary
marks.bi <- array(data = NA, c(n.sps, n.mks, g.len),  # n.sample, n.mark, 0-1 code 
                  dimnames = list(samples, marks, c(1:sum(chrom.leng))))
marks.bi[1, , ] <- rbind(myBinary(cbellum.h3k4me1), myBinary(cbellum.h3k4me3), myBinary(cbellum.h3k27ac), myBinary(cbellum.h3k27me3), myBinary(cbellum.ctcf), myBinary(cbellum.pol2))
marks.bi[2, , ] <- rbind(myBinary(esb4.h3k4me1), myBinary(esb4.h3k4me3), myBinary(esb4.h3k27ac), myBinary(esb4.h3k27me3), myBinary(esb4.ctcf), myBinary(esb4.pol2))
marks.bi[3, , ] <- rbind(myBinary(heart.h3k4me1), myBinary(heart.h3k4me3), myBinary(heart.h3k27ac), myBinary(heart.h3k27me3), myBinary(heart.ctcf), myBinary(heart.pol2))
marks.bi[4, , ] <- rbind(myBinary(kidney.h3k4me1), myBinary(kidney.h3k4me3), myBinary(kidney.h3k27ac), myBinary(kidney.h3k27me3), myBinary(kidney.ctcf), myBinary(kidney.pol2))
marks.bi[5, , ] <- rbind(myBinary(liver.h3k4me1), myBinary(liver.h3k4me3), myBinary(liver.h3k27ac), myBinary(liver.h3k27me3), myBinary(liver.ctcf), myBinary(liver.pol2))
marks.bi[6, , ] <- rbind(myBinary(mel.h3k4me1), myBinary(mel.h3k4me3), myBinary(mel.h3k27ac), myBinary(mel.h3k27me3), myBinary(mel.ctcf), myBinary(mel.pol2))
marks.bi[7, , ] <- rbind(myBinary(smint.h3k4me1), myBinary(smint.h3k4me3), myBinary(smint.h3k27ac), myBinary(smint.h3k27me3), myBinary(smint.ctcf), myBinary(smint.pol2))
marks.bi[8, , ] <- rbind(myBinary(spleen.h3k4me1), myBinary(spleen.h3k4me3), myBinary(spleen.h3k27ac), myBinary(spleen.h3k27me3), myBinary(spleen.ctcf), myBinary(spleen.pol2))
marks.bi[9, , ] <- rbind(myBinary(testis.h3k4me1), myBinary(testis.h3k4me3), myBinary(testis.h3k27ac), myBinary(testis.h3k27me3), myBinary(testis.ctcf), myBinary(testis.pol2))
marks.bi[10, , ] <- rbind(myBinary(thymus.h3k4me1), myBinary(thymus.h3k4me3), myBinary(thymus.h3k27ac), myBinary(thymus.h3k27me3), myBinary(thymus.ctcf), myBinary(thymus.pol2))
marks.bi[11, , ] <- rbind(myBinary(brain.h3k4me1), myBinary(brain.h3k4me3), myBinary(brain.h3k27ac), myBinary(brain.h3k27me3), myBinary(brain.ctcf), myBinary(brain.pol2))
save(marks.bi, file = "~/Dropbox/Network/R2/mark1.rdt")

#--- II. Signal
marks.si <- array(data = NA, c(n.sps, n.mks, g.len),  # n.sample, n.mark, signal 
                  dimnames = list(samples, marks, c(1:sum(chrom.leng))))
marks.si[1, , ] <- rbind(mySignal(cbellum.h3k4me1), mySignal(cbellum.h3k4me3), mySignal(cbellum.h3k27ac), mySignal(cbellum.h3k27me3), mySignal(cbellum.ctcf), mySignal(cbellum.pol2))
marks.si[2, , ] <- rbind(mySignal(esb4.h3k4me1), mySignal(esb4.h3k4me3), mySignal(esb4.h3k27ac), mySignal(esb4.h3k27me3), mySignal(esb4.ctcf), mySignal(esb4.pol2))
marks.si[3, , ] <- rbind(mySignal(heart.h3k4me1), mySignal(heart.h3k4me3), mySignal(heart.h3k27ac), mySignal(heart.h3k27me3), mySignal(heart.ctcf), mySignal(heart.pol2))
marks.si[4, , ] <- rbind(mySignal(kidney.h3k4me1), mySignal(kidney.h3k4me3), mySignal(kidney.h3k27ac), mySignal(kidney.h3k27me3), mySignal(kidney.ctcf), mySignal(kidney.pol2))
marks.si[5, , ] <- rbind(mySignal(liver.h3k4me1), mySignal(liver.h3k4me3), mySignal(liver.h3k27ac), mySignal(liver.h3k27me3), mySignal(liver.ctcf), mySignal(liver.pol2))
marks.si[6, , ] <- rbind(mySignal(mel.h3k4me1), mySignal(mel.h3k4me3), mySignal(mel.h3k27ac), mySignal(mel.h3k27me3), mySignal(mel.ctcf), mySignal(mel.pol2))
marks.si[7, , ] <- rbind(mySignal(smint.h3k4me1), mySignal(smint.h3k4me3), mySignal(smint.h3k27ac), mySignal(smint.h3k27me3), mySignal(smint.ctcf), mySignal(smint.pol2))
marks.si[8, , ] <- rbind(mySignal(spleen.h3k4me1), mySignal(spleen.h3k4me3), mySignal(spleen.h3k27ac), mySignal(spleen.h3k27me3), mySignal(spleen.ctcf), mySignal(spleen.pol2))
marks.si[9, , ] <- rbind(mySignal(testis.h3k4me1), mySignal(testis.h3k4me3), mySignal(testis.h3k27ac), mySignal(testis.h3k27me3), mySignal(testis.ctcf), mySignal(testis.pol2))
marks.si[10, , ] <- rbind(mySignal(thymus.h3k4me1), mySignal(thymus.h3k4me3), mySignal(thymus.h3k27ac), mySignal(thymus.h3k27me3), mySignal(thymus.ctcf), mySignal(thymus.pol2))
marks.si[11, , ] <- rbind(mySignal(brain.h3k4me1), mySignal(brain.h3k4me3), mySignal(brain.h3k27ac), mySignal(brain.h3k27me3), mySignal(brain.ctcf), mySignal(brain.pol2))
save(marks.si, file = "~/Dropbox/Network/R2/mark2.rdt")

#--- Encode: RNA-seq
name1 <- list.files(path = "~/Dropbox/Network/RSEM2",
                    pattern = "*.genes.results")
name2 <- gsub(".genes.results", "", name1)
name3 <- unique(gsub("Rep.*", "", name1))  # sample name
for (i in 1:length(name1)) {
  filepath <- file.path("~/Dropbox/Network/RSEM2", name1[i])
  cat(i, "/", length(name1), name1[i], "\n")
  assign(name2[i], read.delim(filepath)[, -2])
}

rna1 <- matrix(nrow = nrow(get(name2[1])), ncol = length(name2),  # rna: TPM 
               dimnames = list(get(name2[1])$gene_id, name2))
rna2 <- matrix(nrow = nrow(get(name2[1])), ncol = length(name3),  # rna: TPM 
               dimnames = list(get(name2[1])$gene_id, name3))

for (i in 1:length(name2)) {  # TPM with replicates
  rna1[, i] <- get(name2[i])$TPM
}

for (i in 1:length(name3)) {  # TPM mean
  rna2[, i] <- (get(paste(name3[i], "Rep1", sep = ""))$TPM + get(paste(name3[i], "Rep2", sep = ""))$TPM) / 2
}

# heatmap(cor(rna, method = "pearson"))

rna1 <- as.data.frame(rna1)
rna2 <- as.data.frame(rna2)
ens2sym <- read.delim("~/Dropbox/X/ensembl2symbol.map", header = F, row.names = 1)
rna1$symbol <- ens2sym[rownames(rna1), ]
rna2$symbol <- ens2sym[rownames(rna2), ]
rna1 <- aggregate(. ~ symbol, rna1, sum)
rna2 <- aggregate(. ~ symbol, rna2, sum)
rna1 <- data.frame(row.names = rna1$symbol, rna1[, 2:ncol(rna1)])
rna2 <- data.frame(row.names = rna2$symbol, rna2[, 2:ncol(rna2)])

rna3 <- rna1[rowSums(rna1) > ncol(rna1), ]  # FT1: miminum expr
rna3 <- rna3[apply(rna3, 1, function (x) {max(x) > 10}), ]  # FT2: maximal expr

group <- gsub("Rep[12]", "", colnames(rna3))

aov.pval <- rep(NA, nrow(rna3))
for (i in 1:nrow(rna3)) {  # FT3: ANOVA
  if (i %% 2e3 == 0) cat(i, "/", nrow(rna3), "\n")
  dat1 <- data.frame(tpm = as.matrix(rna3)[i, ], group)
  aov1 = aov(tpm ~ group, data = dat1)
  aov.pval[i] <- min(summary(aov1)[[1]][["Pr(>F)"]], na.rm = T)
}

genes <- rownames(rna3[aov.pval < 0.05, ])  # Select genes
rna4 <- rna2[genes, ]  # TPM mean with selected genes

logit <- apply(t(rna4), 2, function (x) {max(x) > 1e2})
rna.f <- rna4[logit, ]  # expression level
rna.f <- log2(rna.f + 1)  # log2 transformation
logit <- apply(t(rna.f), 2, function (x) {max(x) - min(x) > 2})
rna.f <- rna.f[logit, ]  # expression variation

save(rna.f, file = "~/Dropbox/Network/R2/rna.rdt")
