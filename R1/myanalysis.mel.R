# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration: characterizations
# Rev: March 14, 2014
#----------------------------------------------------------------------------------------------
library(ggplot2)
# source("~/Dropbox/Encode/R/encode.mel.R")
load("~/Dropbox/Encode/R/mel.RData")
source("~/Dropbox/Encode/R/myfunction1.R")

#--- Association between marks and genetic features (on binary data) ----------------
marks.pm.bi <- refgene.pm.bi * mel.marks.bi
marks.tx.bi <- refgene.tx.bi * mel.marks.bi
marks.ex.bi <- refgene.ex.bi * mel.marks.bi
marks.it.bi <- marks.tx.bi - marks.ex.bi
marks.ig.bi <- mel.marks.bi - marks.tx.bi
#----------------------------------------------------------------------------------------------
marks.1 <- apply(mel.marks.bi, 2, sum)  # coverage: number of 1
marks.pm.1 <- apply(marks.pm.bi, 2, sum)
marks.tx.1 <- apply(marks.tx.bi, 2, sum)
marks.ex.1 <- apply(marks.ex.bi, 2, sum)
marks.it.1 <- apply(marks.it.bi, 2, sum)
marks.ig.1 <- apply(marks.ig.bi, 2, sum)
#---- Code below takes significant time ------------------------------
# marks.n <- apply(mel.marks.bi, 2, myCount1)
# marks.pm.n <- apply(marks.pm.bi, 2, myCount1)
# marks.tx.n <- apply(marks.tx.bi, 2, myCount1)
# marks.ex.n <- apply(marks.ex.bi, 2, myCount1)
# marks.it.n <- apply(marks.it.bi, 2, myCount1)
# marks.ig.n <- apply(marks.ig.bi, 2, myCount1)
# marks.l <- apply(mel.marks.bi, 2, myLength)
# marks.pm.l <- apply(marks.pm.bi, 2, myLength)
# marks.tx.l <- apply(marks.tx.bi, 2, myLength)
# marks.ex.l <- apply(marks.ex.bi, 2, myLength)
# marks.it.l <- apply(marks.it.bi, 2, myLength)
# marks.ig.l <- apply(marks.ig.bi, 2, myLength)
# save(marks.n, marks.pm.n, marks.tx.n, marks.ex.n, marks.it.n, marks.ig.n, 
#      file = "~/Dropbox/Encode/R/mel1.rdt")
# save(marks.l, marks.pm.l, marks.tx.l, marks.ex.l, marks.it.l, marks.ig.l, 
#      file = "~/Dropbox/Encode/R/mel2.rdt")
load("~/Dropbox/Encode/R/mel1.rdt")
load("~/Dropbox/Encode/R/mel2.rdt")
#----------------------------------------------------------------------------------------------
coverage.1 <- rbind(exon = marks.ex.1,  # exon
                    promoter = marks.pm.1,  # promoter
                    intron = marks.it.1,  # intron
                    intergenic = marks.ig.1)  # inter-genic 
coverage.1 <- coverage.1 / sum(chrom.length) * 100  # in percentage

coverage.2 <- rbind(exon = marks.ex.n,  # exon
                    promoter = marks.pm.n,  # promoter
                    intron = marks.it.n,  # intron
                    intergenic = marks.ig.n)  # inter-genic

#----------------------------------------------------------------------------------------------
bar.dt1 <- data.frame(value = c(as.vector(coverage.1)), 
                      mark = c(apply(as.matrix(colnames(coverage.1)), 1, function (x) {rep(gsub("MEL-", "", x), 4)})),
                      feature = rep(rownames(coverage.1), 7))
bar.dt1$order <- reorder(bar.dt1$mark, bar.dt1$value)

bar.dt2 <- data.frame(value = c(as.vector(coverage.2)) / 10000, 
                      mark = c(apply(as.matrix(colnames(coverage.2)), 1, function (x) {rep(gsub("MEL-", "", x), 4)})),
                      feature = rep(rownames(coverage.2), 7))
bar.dt2$order <- reorder(bar.dt2$mark, bar.dt2$value)
# bar.dt2$order <- reorder(bar.dt2$feature, bar.dt2$value)
pdf("~/Dropbox/Encode/Figures/bar2.pdf", width = 4)
ggplot(bar.dt2, aes(x = order, y = value, fill = feature)) +  # bar.dt1, bar.dt2
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
# xlab("") + ylab("Coverage (Genome %)") +  # bar.dt1
  xlab("") + ylab("Peak number (x 10000)") +  # bar.dt2
  scale_colour_brewer(palette = "Set1") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

n.max = 100000
features = c("genome", "promoter", "gene", "exon", "intron", "intergenic")
hist.dt1 <- data.frame(value = c(rbind(marks.l, marks.pm.l, marks.tx.l, marks.ex.l, marks.it.l, marks.ig.l)), 
                       mark = c(apply(as.matrix(colnames(marks.l)), 1, function (x) {rep(gsub("MEL-", "", x), 6 * n.max)})),
                       feature = rep(c(apply(as.matrix(features), 1, function (x) {rep(x, n.max)})), 7))
hist.dt1 <- hist.dt1[!is.na(hist.dt1$value), ]

pdf(file = "~/Dropbox/Encode/Figures/hist1.pdf", width = 10, height = 8)
ggplot(hist.dt1, aes(x = value, fill = mark)) + 
# geom_density(alpha = .3) + facet_grid(feature ~ .) +
  geom_bar(stat = "density") + facet_grid(feature ~ .) +
# geom_bar(position = "dodge") + facet_grid(feature ~ .) +
  theme_bw() + xlim(c(0, 15)) +
  xlab("Mark length") + ylab("Density") +
  scale_colour_brewer(palette="Set1") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank())
dev.off()
                       
#--- Enrichment ------------------------------------------------------------------
refgene.it.bi <- refgene.tx.bi - refgene.ex.bi
refgene.ig.bi <- rep(1, length(refgene.tx.bi)) - refgene.tx.bi
enr.tx.1 <- (marks.tx.1 / sum(refgene.tx.bi)) / (marks.1 / sum(chrom.length))
enr.ex.1 <- (marks.ex.1 / sum(refgene.ex.bi)) / (marks.1 / sum(chrom.length))
enr.it.1 <- (marks.it.1 / sum(refgene.it.bi)) / (marks.1 / sum(chrom.length))
enr.ig.1 <- (marks.ig.1 / sum(refgene.ig.bi)) / (marks.1 / sum(chrom.length))
enr.tx.n <- (marks.tx.n / sum(refgene.tx.bi)) / (marks.n / sum(chrom.length))
enr.ex.n <- (marks.ex.n / sum(refgene.ex.bi)) / (marks.n / sum(chrom.length))
enr.it.n <- (marks.it.n / sum(refgene.it.bi)) / (marks.n / sum(chrom.length))
enr.ig.n <- (marks.ig.n / sum(refgene.ig.bi)) / (marks.n / sum(chrom.length))
enr <- rbind(enr.tx.1, enr.ex.1, enr.it.1, enr.ig.1, enr.tx.n, enr.ex.n, enr.it.n, enr.ig.n)
rownames(enr) <- c("Gene-C", "Exon-C", "Intron-C", "Intergene-C",
                   "Gene-N", "Exon-N", "Intron-N", "Intergene-N")

bar.dt3 <- data.frame(value = c(enr), 
                      mark = c(apply(as.matrix(colnames(enr)), 1, function (x) {rep(gsub("MEL-", "", x), 8)})),
                      feature = rep(rownames(enr), 7))
bar.dt3$order <- reorder(bar.dt3$feature, bar.dt3$value)

pdf("~/Dropbox/Encode/Figures/bar3.pdf")
ggplot(bar.dt3, aes(x = order, y = value, fill = mark)) +  # bar.dt1, bar.dt2
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  xlab("") + ylab("Enrichment (X)") +
  scale_colour_brewer(palette = "Set1") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()
