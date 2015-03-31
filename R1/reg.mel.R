# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration: stepwise selection of regressors
# Rev: May 7, 2014

library(MASS)
library(ggplot2)
load("~/Dropbox/Encode/R/mel.RData")
source("~/Dropbox/Encode/R/myfunction1.R")

#--- The Pattern -----------------------------------------------------------------------------------
bi <- c(0:1)
pattern <- expand.grid(h3k4me1 = bi, h3k4me3 = bi, h3k9ac = bi, h3k27ac = bi, h3k27me3 = bi,
                       h3k36me3 = bi, h3k79me2 = bi)

tile.dt1 <- data.frame(value = c(as.matrix(pattern)), 
                       idx = rep(paste("P", rownames(pattern), sep = ""), 7), 
                       mark = c(apply(as.matrix(colnames(pattern)), 1, function (x) {rep(x, 128)})))
pdf(file = "~/Dropbox/Encode/Figures/tile1.pdf", height = 2, width = 20)
ggplot(tile.dt1, aes(x = idx, y = mark, fill = value)) + geom_tile() +
  scale_fill_gradient(low="white", high="black") +
  theme_bw() + xlab("") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        legend.position = "none")
dev.off()

#--- The Data --------------------------------------------------------------------------------------
marks.1.tx <- data.frame(row.names = rownames(gene.bed),  # Mark 1s over transcript
                         apply(mel.marks.bi, 2, function(x) {myCount2(gene.bed, x)}))
# pt.marks.1.tx <- as.data.frame(myPattern(marks.1.tx, pattern))
# save(pt.marks.1.tx, file = "~/Dropbox/Encode/R/pt.dt1.rdt")
load("~/Dropbox/Encode/R/pt.dt1.rdt")
pt.marks.1.tx <- pt.marks.1.tx[, -c(which(colSums(pt.marks.1.tx) < 100))]  # 10, 100

pt.sum1 <- colSums(pt.marks.1.tx)
line.dt1 <- data.frame(value = sort(pt.sum1, decreasing = T))
pdf(file = "~/Dropbox/Encode/Figures/line1.pdf", width = 10)
ggplot(line.dt1, aes(x = c(1:dim(line.dt1)[1]), y = value)) + geom_point(size = 5) +
  theme_bw() + xlab("") + ylab("Pattern number") +
  scale_x_continuous(breaks = c(1:dim(line.dt1)[1]), labels = rownames(line.dt1)) +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"))
dev.off()

#--- The Regressions -------------------------------------------------------------------------------
data1 <- cbind(rna = log2(mel.rna$TPM + 1), marks.1.tx)
data1 <- data1[rowSums(data1) > 1, ]  # FT
reg1 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "+"))), data1)  # additive
reg1 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "*"))), data1)  # complete

reg1.step1 <- lm(rna ~ 1, data1) 
reg1.step2 <- stepAIC(reg1.step1, as.formula(paste("~. +", paste(colnames(data1)[-1], collapse = "+"))),
                      direction="forward")  # additive
reg1.step3 <- stepAIC(reg1.step1, as.formula(paste("~. +", paste(colnames(data1)[-1], collapse = "*"))),
                      direction="forward")  # additive & first order interactive

#-- FT: delete Absent patterns
data2 <- cbind(rna = log2(mel.rna$TPM + 1), pt.marks.1.tx)
#-- FT: delete non-expressed, non-marked genes
data2 <- data2[!((data2$V1 == 1) & (data2$rna == 0)), ]  

reg2.step1 <- lm(rna ~ 1, data2)
reg2.step2 <- stepAIC(reg2.step1, as.formula(paste("~. +", paste(colnames(data2)[-1], collapse = "+"))),
                      direction="forward")  # additive
reg2.step3 <- stepAIC(reg2.step1, as.formula(paste(paste("~. + (", paste(colnames(data2)[-1], collapse = "+")), ")^2")),
                      direction="forward")  # additive & first order interaction 
# reg2.step4 <- stepAIC(reg3.step1, as.formula(paste("~. +", paste(colnames(data2)[-1], collapse = "*"))),
#                       direction="forward")  # Error: C stack usage is too close to the limit
idx <- gsub("V", "", names(reg2.step3$coefficients)[-1])

tile.dt2 <- data.frame(value = c(as.matrix(pattern[idx, ])), 
                       idx = rep(paste("P", rownames(pattern[idx, ]), sep = ""), 7), 
                       mark = c(apply(as.matrix(colnames(pattern[idx, ])), 1, function (x) {rep(x, length(idx))})))
pdf(file = "~/Dropbox/Encode/Figures/tile2.pdf", height = 3)
ggplot(tile.dt2, aes(x = idx, y = mark, fill = value)) + geom_tile() +
  scale_fill_gradient(low="white", high="black") +
  theme_bw() + xlab("") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        legend.position = "none")
dev.off()

data3 <- pt.marks.1.tx[, -1]  # Remove non-markded pattern 
reg3.step2 <- list()
for (i in 1:ncol(data3)) {
  print(i)
  reg3.step1 <- lm(as.formula(paste(colnames(data3)[i], "~ 1")), data3)
  reg3.step2[[i]] <- stepAIC(reg3.step1, as.formula(paste(paste("~. + (", paste(colnames(data3)[-i], collapse = "+")), ")^2")),
                        direction="forward")
}

#  Regression between marks (patterns): Understand the results!!! Rebuild the data frame, break the correlations?
#  Transform the data in ways. How to understand this result?
