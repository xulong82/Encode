# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration: characterizations
# Rev: March 14, 2014

library(parallel)
library(pheatmap)
source("~/Dropbox/Encode/R/myfunction1.R")
source("~/Dropbox/Encode/R/encode.R")

#--- Segment ENCODE broadPeak file -----------------------------------------------------------
marks.bi <- cbind(myBinary(brain.h3k4me1), myBinary(brain.h3k4me3), myBinary(brain.h3k27me3), myBinary(brain.h3k27ac),
                  myBinary(spleen.h3k4me1), myBinary(spleen.h3k4me3), myBinary(spleen.h3k27me3), myBinary(spleen.h3k27ac),
                  myBinary(testis.h3k4me1), myBinary(testis.h3k4me3), myBinary(testis.h3k27me3), myBinary(testis.h3k27ac))
colnames(marks.bi) <- c("Brain/h3k4me1", "Brain/h3k4me3", "Brain/h3k27me3", "Brain/h3k27ac",
                        "Spleen/h3k4me1", "Spleen/h3k4me3", "Spleen/h3k27me3", "Spleen/h3k27ac",
                        "Testis/h3k4me1", "Testis/h3k4me3", "Testis/h3k27me3", "Testis/h3k27ac")
rownames(marks.bi) <- 1:nrow(marks.bi)

marks.sig <- cbind(mySeg(brain.h3k4me1), mySeg(brain.h3k4me3), mySeg(brain.h3k27me3), mySeg(brain.h3k27ac),
                   mySeg(spleen.h3k4me1), mySeg(spleen.h3k4me3), mySeg(spleen.h3k27me3), mySeg(spleen.h3k27ac),
                   mySeg(testis.h3k4me1), mySeg(testis.h3k4me3), mySeg(testis.h3k27me3), mySeg(testis.h3k27ac))
colnames(marks.sig) <- colnames(marks.bi)
rownames(marks.sig) <- 1:nrow(marks.sig)

sort.marks <- c(c(1, 5, 9), 1 + c(1, 5, 9), 2 + c(1, 5, 9), 3 + c(1, 5, 9))

refgene.pm.bi <- myBinary(refgene.pm) 
refgene.tx.bi <- myBinary(refgene.tx)
refgene.ex.bi <- myBinary(refgene.ex)

# save(marks.bi, refgene.pm.bi, refgene.tx.bi, refgene.ex.bi, file = "~/Dropbox/Encode/R/bi.RData")
# load("~/Dropbox/Encode/R/bi.RData")
# save(marks.sig, file = "~/Dropbox/Encode/R/marks.sig.RData")
# load("~/Dropbox/Encode/R/marks.sig.RData")

#--- Association between marks and genetic features (on binary data) ----------------
marks.pm.bi <- refgene.pm.bi * marks.bi
marks.tx.bi <- refgene.tx.bi * marks.bi
marks.ex.bi <- refgene.ex.bi * marks.bi

marks.1 <- apply(marks.bi, 2, sum)  # coverage: number of 1
marks.pm.1 <- apply(marks.pm.bi, 2, sum)
marks.tx.1 <- apply(marks.tx.bi, 2, sum)
marks.ex.1 <- apply(marks.ex.bi, 2, sum)

save(marks.bi, marks.pm.bi, marks.tx.bi, marks.ex.bi, file = "~/Dropbox/Encode/R/tohpc.RData")
# souce("~/Dropbox/Encode/R/myhpc.R")  # Number of consecutive 1
load("~/Dropbox/Encode/R/outhpc.RData")

coverage.1 <- rbind(exon = marks.ex.1,  # exon
                    intron = marks.tx.1 - marks.ex.1,  # intron
                    intergenic = marks.1 - marks.tx.1)  # inter-genic 
coverage.1 <- coverage.1 / sum(chrom.length) * 100  # in percentage

coverage.2 <- rbind(exon = marks.ex.nb,  # exon
                    intron = marks.tx.nb - marks.ex.nb,  # intron
                    intergenic = marks.nb - marks.tx.nb)  # inter-genic

par(mfrow = c(1, 1), mar = c(5, 10, 4, 2), font = 2)
barplot(coverage.1[, sort.marks], space = c(.6, .2, .2), horiz = T, las = 1, font = 2,  
        xlab = "Genome %", xlim = c(0, 8), font.lab = 2, col = terrain.colors(3))
barplot(coverage.2[, sort.marks], space = c(.6, .2, .2), horiz = T, las = 1, font = 2,  
        xlab = "Peak Number", xlim = c(0, 140000), font.lab = 2, col = terrain.colors(3))
abline(v = 0, lwd = 1, col = "black")
legend("topright", c("Exon", "Intron", "Non-gene"), fill = colors, box.lwd = 0)

#--- Enrichment ------------------------------------------------------------------
enr.tx.1 <- (marks.tx.1 / sum(refgene.tx.bi)) / (marks.1 / sum(chrom.length))
enr.ex.1 <- (marks.ex.1 / sum(refgene.ex.bi)) / (marks.1 / sum(chrom.length))
enr.tx.nb <- (marks.tx.nb / sum(refgene.tx.bi)) / (marks.nb / sum(chrom.length))
enr.ex.nb <- (marks.ex.nb / sum(refgene.ex.bi)) / (marks.nb / sum(chrom.length))
enr <- rbind(enr.tx.1, enr.ex.1, enr.tx.nb, enr.ex.nb)
rownames(enr) <- c("tx.1", "ex.1", "tx.nb", "ex.nb")

enr.box <- t(enr[, sort.marks])
enr.box <- cbind(enr.box[1:3, ], enr.box[4:6, ], enr.box[7:9, ], enr.box[10:12, ])
enr.box <- enr.box[, c(1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15, 4, 8, 12, 16)]
rownames(enr.box) <- c("Brain", "Spleen", "Testis")
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2), font.lab = 2, font = 2)
boxplot(enr.box, las = 2, xaxt = "n", ylab = "Enrichment (X)", col = terrain.colors(4), font.lab = 2)
legend("topleft", c("h3k4me1", "h3k4me3", "h3k27me3", "h3k27ac"), fill = terrain.colors(4), box.lwd = 0)
axis(1, at = 2.5 + 4 * c(0:3), font = 2, labels = c("TX/coverage", "Exon/coverage", "TX/number", "Exon/number"))

#--- Diversity --------------------------------------------------------------------------
marks.bi.1 <- marks.bi[rowSums(marks.bi) > 0, ]  # Genome with > 1 mark
pcor.marks <- cor(marks.bi.1, method = "pearson")
pheatmap(pcor.marks, fontsize = 16, fontsize_number = 10, display_number = T, 
         treeheight_row = 150, treeheight_col = 0, cellwidth = 38, cellheight = 40)
# dist.marks <- dist(t(marks.bi.1))  # Similar results with Pearson's
# plot(hclust(dist.marks), hang = -1)

#--- Associations with gene expression ---------------------------------------------------
rna.pm.bi <- cbind(myBinary(brain.rna.pm), myBinary(spleen.rna.pm), myBinary(testis.rna.pm))
colnames(rna.pm.bi) <- c("Brain/RNA", "Spleen/RNA", "Testis/RNA")

marks.pm.rna.bi <- cbind(marks.pm.bi[, 1:4] * rna.pm.bi[, 1],  # marked, expressed
  marks.pm.bi[, 5:8] * rna.pm.bi[, 2], marks.pm.bi[, 9:12] * rna.pm.bi[, 3])

marks.pm.rna.cm.bi <- marks.pm.rna.bi[, 1:4] * marks.pm.rna.bi[, 5:8] * marks.pm.rna.bi[, 9:12]
colnames(marks.pm.rna.cm.bi) <- c("h3k4me1", "h3k4me3", "h3k27me3", "h3k27ac")

rna.pm.nb <- apply(rna.pm.bi, 2, myCount1)  # expressed, with > 2 RPKM
marks.pm.rna.nb <- myParallel1(marks.pm.rna.bi)  # marked, expressed
marks.pm.rna.cm.nb <- myParallel1(marks.pm.rna.cm.bi)  # marked, expressed, common

marks.rna.nb <- rbind(rep(marks.pm.rna.cm.nb, 3),  # marked, expressed, common
                      marks.pm.rna.nb - rep(marks.pm.rna.cm.nb, 3),  # marked, expressed, not-common
                      marks.pm.nb - marks.pm.rna.nb,  # marked, not-expressed
                      rna.pm.nb - marks.pm.rna.nb)  # expressed, not marked

ratio1 <- (marks.rna.nb[1, ] / colSums(marks.rna.nb[1:2, ])) * 100

par(mfrow = c(1, 1), mar = c(5, 10, 4, 2), font = 2)
barplot(marks.rna.nb[c(1, 2, 3), sort.marks], space = c(.6, .2, .2), horiz = T, las = 1, font = 2,  
        xlab = "Gene Number", xlim = c(0, 12000), font.lab = 2, col = terrain.colors(3))
abline(v = 0, lwd = 1, col = "black")
legend("topright", c("marked/expressed/common", "marked/expressed/not-common", "marked/not-expressed"), 
       fill = terrain.colors(3), box.lwd = 0)

par(mfrow = c(1, 1), mar = c(10, 20, 4, 20), font = 2)
barplot(ratio1[sort.marks], space = c(.9, .3, .3), las = 2, font = 2, font.lab = 2, 
        ylab = "Ratio (%)", col = terrain.colors(3))
abline(h = 0, lwd = 1, col = "black")

save.image("~/Dropbox/Encode/R/analysis1.20140317.RData")
load("~/Dropbox/Encode/R/analysis1.20140317.RData")