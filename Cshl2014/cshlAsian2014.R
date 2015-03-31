#--- Characterization
library(ggplot2)
library(parallel)

rm(list = ls())
load("~/Dropbox/Network/R1/encodex.rdt")
load("~/Dropbox/Network/R1/esb4Df.rdt")  # ESB4
#----------------------------------------------------------------------------------------------

marks.bi <- marks.bi[, c("k4me1", "k4me3")]
marks.bi <- cbind(marks.bi, overlap = as.numeric(marks.bi[, "k4me1"] & marks.bi[, "k4me3"]))

refgene.tx.marks.bi <- refgene.bi[, "tx"] * marks.bi
refgene.pm2k.marks.bi <- refgene.bi[, "pm2k"] * marks.bi
refgene.ex.marks.bi <- refgene.bi[, "ex"] * marks.bi
refgene.it.marks.bi <- refgene.bi[, "it"] * marks.bi

length.all <- apply(marks.bi, 2, sum)
length.tx <- apply(refgene.tx.marks.bi, 2, sum)
length.pm2k <- apply(refgene.pm2k.marks.bi, 2, sum)
length.ex <- apply(refgene.ex.marks.bi, 2, sum)
length.it <- apply(refgene.it.marks.bi, 2, sum)
length <- rbind(PM2k = length.pm2k, Exon = length.ex, Intron = length.it, Intergene = length.all - length.tx)

length <- length / sum(chrom.leng) * 100  # in percentage

Peak <- function(x) {  # compute number of consecutive 1 (peak) in a binary data
  y <- ifelse(x[1] == 0, 0, 1)  # initiate y
  for (i in 2:length(x)) {
    if (i %% 1e5 == 0) cat(i, "in", length(x), "\n")
    turn <- x[i] - x[i-1]
    y <- y + ifelse(turn == 1, 1, 0)
  }
  return(y)
}

mylist <- list()
mylist[[1]] <- marks.bi
mylist[[2]] <- refgene.tx.marks.bi
mylist[[3]] <- refgene.pm2k.marks.bi
mylist[[4]] <- refgene.ex.marks.bi
mylist[[5]] <- refgene.it.marks.bi

y <- mclapply(mylist, function (x) {apply(x, 2, Peak)}, mc.cores = 5)
y1 <- matrix(unlist(y), nrow = 5, byrow = T)

number <- rbind(PM2k = y[[3]], Exon = y[[4]], Intron = y[[5]], Intergene = y[[1]] - y[[2]])
number <- number * 1e-3  #


bar.dt1 <- data.frame(value = c(as.vector(length)), 
                      mark = rep(colnames(length), each = nrow(length)),
                      feature = rep(rownames(length), ncol(length)))
bar.dt2 <- data.frame(value = c(as.vector(number)), 
                      mark = rep(colnames(number), each = nrow(number)),
                      feature = rep(rownames(number), ncol(number)))

# pdf("~/Dropbox/Network/CSHL2014/char1.pdf", height = 5)
pdf("~/Dropbox/Network/CSHL2014/char2.pdf", height = 5)
# ggplot(bar.dt1, aes(x = mark, y = value, fill = feature)) +  # bar.dt1, bar.dt2
ggplot(bar.dt2, aes(x = mark, y = value, fill = feature)) +  # bar.dt1, bar.dt2
# geom_bar(stat = "identity") +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
# xlab("") + ylab("Genome %") +  # bar.dt1
  xlab("") + ylab("Peak number x 1e-3") +  # bar.dt2
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold", color = "grey30"),
        axis.text.y = element_text(size = 10, face = "bold", color = "grey30"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

#--- enrichment
enr.length.tx <- (length.tx / sum(refgene.bi[, "tx"])) / (length.all / sum(chrom.leng))
enr.length.pm2k <- (length.pm2k / sum(refgene.bi[, "pm2k"])) / (length.all / sum(chrom.leng))
enr.length.ex <- (length.ex / sum(refgene.bi[, "ex"])) / (length.all / sum(chrom.leng))
enr.length.it <- (length.it / sum(refgene.bi[, "it"])) / (length.all / sum(chrom.leng))
enr.length.ig <- ((length.all - length.tx) / sum(1 - refgene.bi[, "tx"])) / (length.all / sum(chrom.leng))

enr.number.pm2k <- (y[[3]] / sum(refgene.bi[, "pm2k"])) / (y[[1]] / sum(chrom.leng))
enr.number.ex <- (y[[4]] / sum(refgene.bi[, "ex"])) / (y[[1]] / sum(chrom.leng))
enr.number.it <- (y[[5]] / sum(refgene.bi[, "it"])) / (y[[1]] / sum(chrom.leng))
enr.number.ig <- ((y[[1]] - y[[2]]) / sum(1 - refgene.bi[, "tx"])) / (y[[1]] / sum(chrom.leng))

enrich.dt <- rbind(PM2k.Length = enr.length.pm2k, 
                   Exon.Length = enr.length.ex,
                   Intron.Length = enr.length.it,
                   Intergene.Length = enr.length.ig,
                   PM2k.Number = enr.number.pm2k,
                   Exon.Number = enr.number.ex,
                   Intron.Number = enr.number.it,
                   Intergene.Number = enr.number.ig)

bar.dt3 <- data.frame(value = c(as.vector(enrich.dt)), 
                      mark = rep(colnames(enrich.dt), each = nrow(enrich.dt)),
                      feature = rep(rownames(enrich.dt), ncol(enrich.dt)))
bar.dt3$type1 <- gsub("\\..*", "", bar.dt3$feature)
bar.dt3$type2 <- gsub("^.*\\.", "", bar.dt3$feature)

pdf("~/Dropbox/Network/CSHL2014/char3.pdf", height = 8)
ggplot(bar.dt3, aes(x = mark, y = value, fill = type1)) +  # bar.dt1, bar.dt2
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(type2 ~ .) +
  theme_bw() + 
  xlab("") + ylab("Enrichment") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold", color = "grey30"),
        axis.text.y = element_text(size = 10, face = "bold", color = "grey30"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

#--- correlation & association
data0 <- esb4.df
data1 <- data0[apply(data0[, grep("length", colnames(data0))], 1, function (x) {sum(as.logical(x)) > 30}), ]

# rna1 <- data1$Esb4
rna2 <- scale(data1$Esb4)  # Z

geno <- log2(as.matrix(data1[, -c(1:11)]) + 1)  # log2
geno.sd <- apply(geno, 2, sd)
geno <- geno[, geno.sd > 0]
geno <- geno[, grep("k4me", colnames(geno))]
geno <- geno[, -grep("u[35]\\.", colnames(geno))]
geno <- geno[, -grep("en\\.", colnames(geno))]
geno <- geno[, -grep("pm5k\\.", colnames(geno))]
geno <- geno[, -grep("tss1k\\.", colnames(geno))]
geno <- geno[, -grep("\\.peak", colnames(geno))]
geno <- geno[, -grep("\\.signal", colnames(geno))]
geno <- cbind(geno[, grep("intensity", colnames(geno))], 
              k4me1.pm2k.intensity = geno[, "k4me1.pm2k.length"] / 10,
              k4me3.pm2k.intensity = geno[, "k4me3.pm2k.length"] / 10)

cat(date(), "correlation on single epigenotype \n")
single <- matrix(nrow = ncol(geno), ncol = 5, 
                 dimnames = list(colnames(geno), c("estimate", "p.value", "beta", "logLik", "r2")))
for (i in 1:ncol(geno)) {  # per epigenotype
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

cat(date(), "linear regression for epigenotype pair \n")
pairs <- combn(colnames(geno), 2)
pair <- matrix(0, nrow = ncol(pairs), ncol = 5, 
  dimnames = list(paste("pair", 1:ncol(pairs), sep = ""), c("interaction", "se", "p.value", "r2", "logLik")))
for (i in 1: ncol(pairs)) {
  if (i %% 1e3 == 0) cat(i, "in", ncol(pairs), "\n")
# data1 <- data.frame(rna = rna1, geno[, pairs[, i]])
  data1 <- data.frame(rna = rna2, geno[, pairs[, i]])
  model1 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "+"))), data1)  # additive
  model2 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "*"))), data1)  # complete
  if (nrow(summary(model2)$coefficients) == 4) 
    pair[i, ] <- c(summary(model2)$coefficients[, "Estimate"][4], 
                   summary(model2)$coefficients[, "Std. Error"][4],
                   summary(model2)$coefficients[, "Pr(>|t|)"][4], 
                   summary(model2)$r.squared, logLik(model2) - logLik(model1))
}

feature.map <- data.frame(row.names = c("tx", "en", "tss1k", "pm2k", "pm5k", "ex", "it", "u5", "u3"),
                          label = c("Transcript", "Enhancer", "TSS (1K)", "Promoter (2K)", "Promoter (5K)", "Exon", "Intron", "5'UTR", "3'UTR"))

#--- effect / pvalue
single.dt <- as.data.frame(single)
single.dt$mark <- gsub("\\..*", "", rownames(single.dt))
single.dt$mark <- factor(single.dt$mark, levels = unique(single.dt$mark))
single.dt$feature <- gsub("^.*\\.(.*?)", "", rownames(single.dt))
single.dt$feature <- gsub("\\..*", "", single.dt$feature)
single.dt$feature <- feature.map[single.dt$feature, "label"]
single.dt$q.value <- p.adjust(single.dt$p.value, method = "fdr")
single.dt$FDR <- ifelse(single.dt$q.value < .05, "< 0.05", "> 0.05")
single.dt$FDR <- factor(single.dt$FDR, levels = c("> 0.05", "< 0.05"))
pdf("~/Dropbox/Network/CSHL2014/singleCor.pdf", height = 5)
ggplot(single.dt, aes(x = mark, y = estimate, color = feature)) +  
  geom_point(aes(shape = FDR), size = 5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_bw() + xlab("") + ylab("Effect") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, color = "gray30", face = "bold"),
        axis.text.y = element_text(size = 10, color = "gray30", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) + 
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold", color = "grey30"), legend.key = element_blank()) 
dev.off()

pair.dt <- as.data.frame(pair)
pair.dt$q.value <- p.adjust(pair.dt$p.value, method = "fdr")
pair.dt$significance <- -log10(pair.dt$p.value)

pair.dt$label <- apply(pairs, 2, function (x) {paste(x, collapse = " X ")})
pair.dt$label <- gsub("\\.intensity", "", pair.dt$label)
pair.dt$FDR <- ifelse(pair.dt$q.value < .05, "< 0.05", "> 0.05")
pair.dt$FDR <- factor(pair.dt$FDR, levels = c("> 0.05", "< 0.05"))

limits <- aes(ymax = resp + se, ymin=resp - se)

pdf("~/Dropbox/Network/CSHL2014/pairLm.pdf", height = 8, width = 10)
ggplot(pair.dt, aes(x = label, y = interaction)) +  
  geom_point(aes(shape = FDR, color = FDR), size = 6) + 
  geom_errorbar(aes(ymax = interaction + se, ymin = interaction - se), width = 0.2) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_bw() + xlab("") + ylab("Effect") +
  scale_colour_manual(values = c("darkred", "darkgreen")) +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = -90, color = "grey30"),
        axis.text.y = element_text(size = 10, face = "bold", color = "grey30"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold", color = 'grey30'), legend.key = element_blank()) 
dev.off()

#--- INTERACTION PLOT
features <- c("k4me1.tx.intensity", "k4me3.tx.intensity")

data3 <- data.frame(rna = rna2, geno[, features])
data3$f1 <- data3$f2 <- rep("q2", nrow(data3))
data3$f1[data3[, 2] < quantile(data3[, 2], .25)] <- "q1"
data3$f1[data3[, 2] > quantile(data3[, 2], .75)] <- "q3"
data3$f2[data3[, 3] < quantile(data3[, 3], .25)] <- "q1"
data3$f2[data3[, 3] > quantile(data3[, 3], .75)] <- "q3"

data3$group <- paste(data3$f1, data3$f2, sep = "/")
data4 <- data3[-grep("q2", data3$group), ]
table(data4$group)

data4$group <- gsub("q1/q1", "me1-/m3-", data4$group)
data4$group <- gsub("q1/q3", "me1-/m3+", data4$group)
data4$group <- gsub("q3/q1", "me1+/m3-", data4$group)
data4$group <- gsub("q3/q3", "me1+/m3+", data4$group)
pdf("~/Dropbox/Network/CSHL2014/q1q3plot1.pdf", width = 6, height = 8)
ggplot(data4, aes(x = group, y = rna, fill = group)) + geom_boxplot() +
# geom_hline(yintercept = 65, color = "darkred", size = 2, linetype = 2) + 
  theme_bw() + xlab("") + ylab("RNA") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 0, face = "bold", angle = -90, color = "grey30"),
        axis.text.y = element_text(size = 10, face = "bold", color = "grey30"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

model1 <- lm(as.formula(paste("rna ~ ", paste(colnames(data3)[2], collapse= "+"))), data3)  # additive
model2 <- lm(as.formula(paste("rna ~ ", paste(colnames(data3)[3], collapse= "+"))), data3)  # additive
model3 <- lm(as.formula(paste("rna ~ ", paste(colnames(data3)[2:3], collapse= "+"))), data3)  # additive
model4 <- lm(as.formula(paste("rna ~ ", paste(colnames(data3)[2:3], collapse= "*"))), data3)  # complete

model1.dt <- data.frame(feature = names(model1$coefficients), effect = model1$coefficients, se = summary(model1)$coefficients[, "Std. Error"])
model1.dt$model <- rep("model1", nrow(model1.dt))
model2.dt <- data.frame(feature = names(model2$coefficients), effect = model2$coefficients, se = summary(model2)$coefficients[, "Std. Error"])
model2.dt$model <- rep("model2", nrow(model2.dt))
model3.dt <- data.frame(feature = names(model3$coefficients), effect = model3$coefficients, se = summary(model3)$coefficients[, "Std. Error"])
model3.dt$model <- rep("model3", nrow(model3.dt))
model4.dt <- data.frame(feature = names(model4$coefficients), effect = model4$coefficients, se = summary(model4)$coefficients[, "Std. Error"])
model4.dt$model <- rep("model4", nrow(model4.dt))
model.dt <- rbind(model1.dt, model2.dt, model3.dt, model4.dt)

model.dt$feature <- gsub("\\(Intercept\\)", "Intercept", model.dt$feature)
model.dt$position <- as.numeric(gsub("model", "", model.dt$model))
model.dt$feature <- factor(model.dt$feature, levels = unique(model.dt$feature))
model.dt$feature <- gsub("\\.intensity", "", model.dt$feature)

pdf("~/Dropbox/Network/CSHL2014/model1.pdf", width = 5)
ggplot(model.dt, aes(x = position, y = effect, color = feature)) +  
  geom_point(size = 5) + geom_line(linetype = 3, width = 2) +
  geom_errorbar(aes(ymax = effect + se, ymin = effect - se), width = 0.2) + 
  theme_bw() + xlab("") + ylab("Effect") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()
