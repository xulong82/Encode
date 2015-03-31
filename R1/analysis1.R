library(MASS)

rm(list = ls())
load("~/Dropbox/Network/R1/esb4Df.rdt")  # ESB4
load("~/Dropbox/Network/R1/melDf.rdt")  # MEL
load("~/Dropbox/Network/R1/liverDf.rdt")  # Liver

data0 <- esb4.df
data0 <- mel.df
data0 <- liver.df

#--- selection
# data1 <- esb4.df[esb4.df$rna > 5, ]  # choose on rna
data1 <- data0[apply(data0[, grep("length", colnames(data0))], 1, function (x) {sum(as.logical(x)) > 30}), ]

# rna1 <- data1$Esb4
rna2 <- scale(data1$Esb4)  # Z
rna2 <- scale(data1$Mel)  # Z
rna2 <- scale(data1$Liver)  # Z

geno <- log2(as.matrix(data1[, -c(1:11)]) + 1)  # log2: CAVEAT: DO I LOG2 ALL FEATURES, OR ANY THE LENGTH
geno.sd <- apply(geno, 2, sd)
geno <- geno[, geno.sd > 0]

#--- correlation & association
cat(date(), "correlation on single epigenotype \n")
single <- matrix(nrow = ncol(geno), ncol = 5, 
  dimnames = list(colnames(geno), c("estimate", "p.value", "beta", "logLik", "r2")))
for (i in 1:ncol(geno)) {  # per epigenotype
  if (i %% 1e2 == 0) cat(i, "in", ncol(geno), "\n")
  model1 <- lm(rna2 ~ 1)
  model2 <- lm(rna2 ~ geno[, i])
  cor1 = cor.test(rna2, geno[, i], method = "pearson")
# cor1 = cor.test(rna2, geno[, i], method = "kendall")
  single[i, "estimate"] = cor1$estimate
  single[i, "p.value"] = cor1$p.value
  single[i, "beta"]  = model2$coefficients[2]
  single[i, "logLik"] = logLik(model2) - logLik(model1)
  single[i, "r2"] = summary(model2)$r.squared
}

sort(single[, "estimate"], decreasing = T)
sort(abs(single[, "estimate"]), decreasing = T)
sort(abs(single[, "p.value"]))
sort(single[, "r2"], decreasing = T)

cat(date(), "linear regression for epigenotype pair \n")
pairs <- combn(colnames(geno), 2)
pair <- matrix(0, nrow = ncol(pairs), ncol = 7, 
  dimnames = list(paste("pair", 1:ncol(pairs), sep = ""), 
  c("intercept", "mark1", "mark2", "interaction", "p.value", "r2", "logLik")))
for (i in 1: ncol(pairs)) {
  if (i %% 1e3 == 0) cat(i, "in", ncol(pairs), "\n")
  # data1 <- data.frame(rna = rna1, geno[, pairs[, i]])
  data1 <- data.frame(rna = rna2, geno[, pairs[, i]])
  model1 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "+"))), data1)  # additive
  model2 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "*"))), data1)  # complete
  if (nrow(summary(model2)$coefficients) == 4) 
    pair[i, ] <- c(summary(model2)$coefficients[, "Estimate"], 
                   summary(model2)$coefficients[, "Pr(>|t|)"][4], 
                   summary(model2)$r.squared, logLik(model2) - logLik(model1))
  # anova(model)
}

summary(pair[, "intercept"])
summary(pair[, "mark1"])
summary(pair[, "mark2"])
summary(pair[, "interaction"])
summary(pair[, "r2"])
summary(pair[, "logLik"])
hist(pair[, "r2"])

#---
cat(date(), " run stepwise linear regression \n")
aic.step1 <- lm(rna2 ~ 1, data3) 
aic.step2 <- stepAIC(aic.step1, as.formula(paste("~. +", paste(colnames(geno), collapse = "+"))),
                     direction="forward")  # additive
aic.step2 <- stepAIC(aic.step1, as.formula(paste(paste("~. + (", paste(colnames(geno), collapse = "+")), ")^2")),
                     direction="forward")  # additive & first order interaction 
aic.step2 <- stepAIC(aic.step1, as.formula(paste("~. +", paste(colnames(geno), collapse = "*"))),
                     direction="forward")  # additive & first order interactive
mod.aic <- aic.step2

summary(mod.aic)

#---
save(single, pair, pairs, file = "~/Dropbox/Network/R1/esb4Analysis.rdt")
save(single, pair, pairs, file = "~/Dropbox/Network/R1/melAnalysis.rdt")
save(single, pair, pairs, file = "~/Dropbox/Network/R1/liverAnalysis.rdt")
