# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - linear regression
# Rev: March 13, 2014

library(pheatmap)
source("~/Dropbox/Encode/R/myfunction2.R")

load("~/Dropbox/Encode/R/bi.RData")  # marks.bi, refgene.[ex|pm|tx].bi
load("~/Dropbox/Encode/R/marks.sig.RData")  # marks.sig
load("~/Dropbox/Encode/R/rna.RData")  # [brain|spleen|testis].rna, RPKM > 2

#---------------------------------------------------------------------------------------------------
gene.space <- unique(rbind(brain.rna[, 1:4], spleen.rna[, 1:4], testis.rna[, 1:4]))
gene.space <- gene.space[!duplicated(gene.space$gene), ]  # Remove duplicated gene items
gene.space <- gene.space[!is.na(gene.space$gene), ]  # Remove NA

brain.rna <- merge(gene.space, aggregate(RPKM ~ gene, data = brain.rna, max), by = "gene", all.x = T)
spleen.rna <- merge(gene.space, aggregate(RPKM ~ gene, data = spleen.rna, max), by = "gene", all.x = T)
testis.rna <- merge(gene.space, aggregate(RPKM ~ gene, data = testis.rna, max), by = "gene", all.x = T)

rna.rpkm <- data.frame(row.names = brain.rna$gene, brain.rna[, 2:4],
                       brain.rna$RPKM, spleen.rna$RPKM, testis.rna$RPKM) 
rna.rpkm[is.na(rna.rpkm)] <- 0

hist(rna.rpkm[, 4], xlim = c(-10, 30), breaks = 10000)
hist(rna.rpkm[rowSums(rna.rpkm[, 4:6]) > 10, 4], xlim = c(-10, 30), breaks = 10000)

gene.info.tx <- rna.rpkm[, 1:3]
gene.info.pm <- data.frame(row.names = rownames(gene.info.tx),
                           chrom = gene.info.tx$chrom, 
                           chromStart = gene.info.tx$chromStart - 2000, 
                           chromEnd = gene.info.tx$chromStart + 2000) 

marks.ex.bi <- refgene.ex.bi * marks.bi
marks.tx <- data.frame(row.names = rownames(gene.info.tx),  # Mark number over transcript
                       apply(marks.bi, 2, function(x) {myCount2(gene.info.tx, x)}))
marks.pm <- data.frame(row.names = rownames(gene.info.pm),  # Mark number over promoter
                       apply(marks.bi, 2, function(x) {myCount2(gene.info.pm, x)}))
marks.ex <- data.frame(row.names = rownames(gene.info.tx),  # Mark number over exon
                       apply(marks.ex.bi, 2, function(x) {myCount2(gene.info.tx, x)}))
rna.marks <- cbind(rna.rpkm, marks.tx, marks.pm, marks.ex)

marks.ex.sig <- refgene.ex.bi * marks.sig
marks.sig.tx <- data.frame(row.names = rownames(gene.info.tx),  # Mark number over transcript
                           apply(marks.sig, 2, function(x) {myCount2(gene.info.tx, x)}))
marks.sig.pm <- data.frame(row.names = rownames(gene.info.pm),  # Mark number over promoter
                           apply(marks.sig, 2, function(x) {myCount2(gene.info.pm, x)}))
marks.sig.ex <- data.frame(row.names = rownames(gene.info.tx),  # Mark number over exon
                           apply(marks.ex.sig, 2, function(x) {myCount2(gene.info.tx, x)}))
rna.marks.sig <- cbind(rna.rpkm, marks.sig.tx, marks.sig.pm, marks.sig.ex)

save(rna.marks, rna.marks.sig, file = "~/Dropbox/Encode/R/rna.marks.RData")
load("~/Dropbox/Encode/R/rna.marks.RData")

#---------------------------------------------------------------------------------------------------
rna.marks[, 4:6] <- log2(rna.marks[, 4:6] + 1)  # log2 transformation of rpkm
testis.tx <- rna.marks[rowSums(rna.marks[15:18]) > 0, c(6, 15:18)]  # marked
testis.pm <- rna.marks[rowSums(rna.marks[27:30]) > 0, c(6, 27:30)]
testis.ex <- rna.marks[rowSums(rna.marks[39:42]) > 0, c(6, 39:42)]
testis.tx <- testis.tx[testis.tx[, 1] > 0, ]  # expressed
testis.pm <- testis.pm[testis.pm[, 1] > 0, ]
testis.ex <- testis.ex[testis.ex[, 1] > 0, ]
pairs(testis.tx, panel = panel.smooth, main = "TX:B")
pairs(testis.pm, panel = panel.smooth, main = "PM:B")
pairs(testis.ex, panel = panel.smooth, main = "EX:B")

rna.marks.sig[, 4:6] <- log2(rna.marks.sig[, 4:6] + 1)  # log2 transformation of rpkm
testis.sig.tx <- rna.marks.sig[rowSums(rna.marks.sig[15:18]) > 0, c(6, 15:18)]  # marked
testis.sig.pm <- rna.marks.sig[rowSums(rna.marks.sig[27:30]) > 0, c(6, 27:30)]
testis.sig.ex <- rna.marks.sig[rowSums(rna.marks.sig[39:42]) > 0, c(6, 39:42)]
testis.sig.tx <- testis.sig.tx[testis.sig.tx[, 1] > 0, ]  # expressed
testis.sig.pm <- testis.sig.pm[testis.sig.pm[, 1] > 0, ]
testis.sig.ex <- testis.sig.ex[testis.sig.ex[, 1] > 0, ]
pairs(testis.sig.tx, panel = panel.smooth, main = "TX:S")
pairs(testis.sig.pm, panel = panel.smooth, main = "PM:S")
pairs(testis.sig.ex, panel = panel.smooth, main = "EX:S")

summary(testis.tx)
summary(testis.pm)
summary(testis.ex)

myThr <- function(x) {  # Thresholize the data
  x[x < quantile(x, 0.3)] <- 0
  x[x >= quantile(x, 0.3) & x < quantile(x, 0.6)] <- 1
  x[x >= quantile(x, 0.6)] <- 2
  return(x)
}

testis.pm.thr <- as.data.frame(cbind(RPKM = testis.pm[, 1], apply(testis.pm[, 2:5], 2, myThr)))
testis.ex.thr <- as.data.frame(cbind(RPKM = testis.ex[, 1], apply(testis.ex[, 2:5], 2, myThr)))

pairs(testis.pm.thr, panel = panel.smooth, main = "TX:T")

pheatmap(cor(testis.pm), fontsize = 16, fontsize_number = 10, display_number = T, 
         treeheight_row = 150, treeheight_col = 0, cellwidth = 38, cellheight = 40)

model1 <- lm(testis.rna.RPKM ~ 
             Testis.h3k4me1 + Testis.h3k4me3 + Testis.h3k27me3 + Testis.h3k27ac, 
             data = testis.pm)
summary(model1)

model2 <- lm(testis.rna.RPKM ~ 
             h3k4me1 + h3k4me3 + h3k27me3 + h3k27ac + h3k4me1 * h3k27ac + h3k4me3 * h3k27ac, 
             data = testis.lm3)
summary(model2)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 6), font = 2)
mypvalue1 <- coef(summary(model1))[, 4]
mypvalue2 <- coef(summary(model3))[, 4]
plot(-log10(mypvalue1[-1]), las = 1, bty = "l", col = "red", pch = 19,
     xaxt = "n", xlab = "Marks", ylab = "-Log10(p)", ylim = c(0, 50))
     axis(1, at=1:4, labels = c("H3k4me1", "H3k4me3", "H3k27me3", "H3k27ac"))
plot(-log10(mypvalue2[-1]), type = "b", las = 1, bty = "l", col = "red", pch = 19,
     xaxt = "n", xlab = "Marks", ylab = "-Log10(p)", ylim = c(0, 50))
     axis(1, at=1:6, labels = c("H3k4me1", "H3k4me3", "H3k27me3", "H3k27ac", "\nH3k4me1\nH3k27ac", "\nH3k4me3\nH3k27ac"))
