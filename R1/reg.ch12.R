# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration: stepwise selection of regressors
# Rev: May 14, 2014

library(MASS)
library(ggplot2)
load("~/Dropbox/Encode/R/ch12.RData")
source("~/Dropbox/Encode/R/myfunction1.R")

#--- The Data --------------------------------------------------------------------------------------
marks.1.tx <- data.frame(row.names = rownames(gene.bed),  # Mark 1s over transcript
                         apply(ch12.marks.bi, 2, function(x) {myCount2(gene.bed, x)}))

#--- The Regressions -------------------------------------------------------------------------------
#--- RNA ~ marks -----------------------------------------------------
data1 <- log2(cbind(rna = ch12.rna$TPM, marks.1.tx) + 1)  # Log2 transform
data1 <- data1[!((data1[, 1] < 5) & (rowSums(data1[, -1]) < 10)), ]  # Low expression & Low marking

reg0 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "+"))), data1)  # additive
reg0 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "*"))), data1)  # complete

n.start <- 2
n.end <- 6
n.combo <- (n.end - 1) * (n.end - 2) /2
reg1 <- matrix(nrow = n.combo, ncol = 4, dimnames = list(1:n.combo, c("m1", "m2", "lrt", "pval")))
k <- 0
for (i in n.start:(n.end - 1)) {
  for (j in (i + 1): n.end) {
    k = k + 1
    cat(k, paste(colnames(data1)[i], colnames(data1)[j], sep = ":"), "\n")
    dt.tmp <- data1[, c(1, i, j)]
    g0 <- lm(as.formula(paste("rna ~ ", paste(colnames(dt.tmp)[-1], collapse = "+"))), data = dt.tmp)
    g1 <- lm(as.formula(paste(paste("rna ~ (", paste(colnames(dt.tmp)[-1], collapse = "+")), ")^2")), data = dt.tmp)
#   g0 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse = "+"))), data = data1)
#   g1 <- lm(as.formula(paste(paste("rna ~ ", paste(colnames(data1)[-1], collapse = "+")),
#                             paste("+", paste(colnames(dt.tmp)[-1], collapse = "*")))), data = data1)
    lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
    pval <- 1 - pchisq(lrt, df = 1)  # Return p value
    reg1[k, ] <- c(i, j, lrt, pval)
  }
}  #  LRT extreme large and p extreme low: R2 very low: eventhough model improves significantly (high LRT, 
   #  because the degre of freedom is huge), they are both bad models (R2, larger than .5 is reasonable)

#--- RNA ~ marks: stepAIC --------------------------------------------
reg2.step1 <- lm(rna ~ 1, data1) 
reg2.step2 <- stepAIC(reg2.step1, as.formula(paste("~. +", paste(colnames(data1)[-1], collapse = "+"))),
                      direction="forward")  # additive
reg2.step3 <- stepAIC(reg2.step1, as.formula(paste("~. +", paste(colnames(data1)[-1], collapse = "*"))),
                      direction="forward")  # additive & first order interactive
#--- marks ~ marks ---------------------------------------------------
data2 <- log2(marks.1.tx[rowSums(marks.1.tx) > 128, ] + 1)
n.start <- 2
n.end <- ncol(data2) 
n.combo <- (n.end - 1) * (n.end - 2) /2
reg2 <- array(data = NA, c(n.end, n.combo, 4), 
              dimnames = list(colnames(data2), 1:n.combo, c("m1", "m2", "lrt", "pval")))
for (l in 1:ncol(data2)) {
  k <- 0
  dt.tmp1 <- cbind(y = data2[, l], data2[, -l])
  for (i in n.start:(n.end - 1)) {
    for (j in (i + 1): n.end) {
      k = k + 1
      cat(k, paste(colnames(dt.tmp1)[i], colnames(dt.tmp1)[j], sep = ":"), "\n")
      dt.tmp2 <- dt.tmp1[, c(1, i, j)]
      g0 <- lm(as.formula(paste("y ~ ", paste(colnames(dt.tmp2)[-1], collapse = "+"))), data = dt.tmp2)
      g1 <- lm(as.formula(paste(paste("y ~ (", paste(colnames(dt.tmp2)[-1], collapse = "+")), ")^2")), data = dt.tmp2)
      # g0 <- lm(as.formula(paste("y ~ ", paste(colnames(data1)[-1], collapse = "+"))), data = data1)
      # g1 <- lm(as.formula(paste(paste("y ~ ", paste(colnames(data1)[-1], collapse = "+")),
      #                           paste("+", paste(colnames(dt.tmp)[-1], collapse = "*")))), data = data1)
      lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
      pval <- 1 - pchisq(lrt, df = 1)  # Return p value
      reg2[l, k, ] <- c(i, j, lrt, pval)
    }
  }
}  #  LRT extreme large and p extreme low
#*** A PROBLEM OF LINEAR REGRESSION IN TERMS OF DATA SIZE AND VARIABLE SIZE ***
#--- mark-mark interaction needs all genomic scale, not restricted to the genomic regions
#--- start from genomic regions, map to genetic region is a must do option
data3 <- ch12.marks.bi[rowSums(ch12.marks.bi) >= 2, ]

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
