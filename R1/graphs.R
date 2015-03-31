rm(list = ls())
load("~/Dropbox/Network/R1/esb4Analysis.rdt")
load("~/Dropbox/Network/R1/liverAnalysis.rdt")

#--- scatterplot
pdf("~/Dropbox/Network/R1.fig/pairs.pdf")
pairs(cbind(rna2, geno[, grep("length", colnames(geno))]))
pairs(cbind(rna2, geno[, grep("peak", colnames(geno))]))
pairs(cbind(rna2, geno[, grep("intensity", colnames(geno))]))
pairs(cbind(rna2, geno[, grep("signal", colnames(geno))]))
dev.off()

feature.map <- data.frame(row.names = c("tx", "en", "tss1k", "pm2k", "pm5k", "ex", "it", "u5", "u3"),
  label = c("Transcript", "Enhancer", "TSS (1K)", "Promoter (2K)", "Promoter (5K)", "Exon", "Intron", "5'UTR", "3'UTR"))

#--- effect / pvalue
single.dt <- as.data.frame(single)
# rownames(single.dt) <- gsub("active\\.en", "active", rownames(single.dt))
# rownames(single.dt) <- gsub("inactive\\.en", "inactive", rownames(single.dt))
single.dt$mark <- gsub("\\..*", "", rownames(single.dt))
single.dt$mark <- factor(single.dt$mark, levels = unique(single.dt$mark))
single.dt$type <- gsub("^.*\\.", "", rownames(single.dt))
single.dt$feature <- gsub("^.*\\.(.*?)", "", rownames(single.dt))
single.dt$feature <- gsub("\\..*", "", single.dt$feature)
single.dt$feature <- feature.map[single.dt$feature, "label"]
single.dt$significance <- -log10(single.dt$p.value)
single.dt$q.value <- p.adjust(single.dt$p.value, method = "fdr")
single.dt$FDR <- ifelse(single.dt$q.value < .05, "< 0.05", "> 0.05")
single.dt$FDR <- factor(single.dt$FDR, levels = c("> 0.05", "< 0.05"))
# pdf("~/Dropbox/Network/R1.fig/singleLm.pdf", width = 10)
pdf("~/Dropbox/Network/R1.fig/singleCor.pdf", width = 10)
ggplot(single.dt, aes(x = mark, y = estimate, color = feature)) +  
# ggplot(single.dt, aes(x = mark, y = beta, color = feature)) +  
# geom_point(aes(size = logLik)) + 
  geom_point(aes(shape = FDR)) + 
# geom_point(aes(size = significance, shape = FDR)) + 
  facet_grid(type ~ .) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
# scale_y_continuous(lim = c(-.2, .2)) +
  theme_bw() + xlab("") + ylab("Effect") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, color = "gray30", face = "bold"),
        axis.text.y = element_text(size = 10, color = "gray30", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) + 
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
#       legend.title = element_blank(), legend.key = element_blank()) 
        legend.title = element_text(size = 10, face = "bold", color = "grey30"), legend.key = element_blank()) 
dev.off()


pair.dt <- as.data.frame(pair)
quantile(abs(pair.dt$logLik), .95)
pair.dt$cut <- ifelse(abs(pair.dt$logLik) > quantile(abs(pair.dt$logLik), .95), "95 quantile", "95 quantile")
pair.dt$q.value <- p.adjust(pair.dt$p.value, method = "fdr")
pair.dt$significance <- -log10(pair.dt$p.value)

pair.dt$label <- apply(pairs, 2, function (x) {paste(x, collapse = " X ")})
pair.dt$TF <- rep("NO", nrow(pair.dt))
pair.dt$TF[grep("ctcf", pair.dt$label)] <- "YES"
pair.dt$TF[grep("pol2", pair.dt$label)] <- "YES"
pair.dt$TF <- factor(pair.dt$TF, levels = c("YES", "NO"))
pair.dt$mark1 <- gsub("\\..*", "", pair.dt$label)
pair.dt$mark2 <- gsub("^.*X ", "", pair.dt$label)
pair.dt$mark2 <- gsub("\\..*", "", pair.dt$mark2)
x = t(apply(pair.dt[, 2:3], 1, sort))
pair.dt$pair <- paste(x[, 1], x[, 2], sep = " X ")

pair.dt1 <- pair.dt[pair.dt$q.value < .05, ]
pair.dt1 <- pair.dt1[abs(pair.dt1$interaction) > 1, ]
pair.dt1 <- pair.dt1[pair.dt1$mark1 != pair.dt1$mark2, ]

pdf("~/Dropbox/Network/R1.fig/groups.pdf", width = 12)
par(mar = c(20, 4, 4, 2))
barplot(table(pair.dt1$pair), beside = T, lwd = 3, las = 3, col = "darkgreen", font.axis = 2)
abline(0, 0, lwd = 3, col = "black")
dev.off()

pair.dt2 <- pair.dt1[order(pair.dt1$pair), ]
pair.dt2$order <- 1:nrow(pair.dt2)

group.lab <- pair.dt2$pair
group.lab <- group.lab[!duplicated(group.lab)]
group.len <- rep(0, length(group.lab))
for (i in 1:length(group.lab)) group.len[i] = length(which(pair.dt2$pair == group.lab[i]))
group.pos <- c(0, cumsum(group.len))[-(length(group.len) + 1)] + group.len / 2
group.pos[group.lab == "k36me3 X k4me1"] = group.pos[group.lab == "k36me3 X k4me1"] + 2

pdf("~/Dropbox/Network/R1.fig/pairLm1.pdf", width = 12)
# ggplot(pair.dt2, aes(x = label, y = interaction)) +  
ggplot(pair.dt2, aes(x = order, y = interaction)) +  
  geom_point(aes(color = TF, size = significance)) + 
# geom_point(aes(color = TF, size = logLik)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
# geom_vline(xintercept = which(pair.dt$pair == "k36me3|k9me3"), linetype = "solid", size = 1) +
# geom_point(aes(size = logLik)) + 
  scale_x_continuous(breaks = group.pos, labels = group.lab) +
# scale_y_continuous(lim = c(-100, 80)) +
  theme_bw() + xlab("") + ylab("Effect") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 7, face = "bold", angle = -90, color = "grey30"),
# theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold", color = "grey30"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold", color = 'grey30'), legend.key = element_blank()) 
dev.off()

pdf("~/Dropbox/Network/R1.fig/singleR2.pdf", width = 5)
# pdf("~/Dropbox/Network/R1.fig/pairR2.pdf", width = 5)
ggplot(single.dt, aes(x = r2)) + 
# ggplot(pair.dt, aes(x = r2)) + 
  # ggplot(pair.dt, aes(x = interaction)) + 
  geom_histogram(colour = "white", fill = "darkgreen", binwidth = .01) +
  # geom_histogram(colour = "white", aes(fill = cut), binwidth = .1) +
  scale_fill_manual(values=c("darkgreen", "darkred")) +
  # scale_x_continuous(lim = c(-5, 5)) +
  theme_bw() + xlab("R2") + ylab("Count") +
  # theme_bw() + xlab("Effect") + ylab("Count") +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 18, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

aic.dt <- as.data.frame(summary(mod.aic)$coefficients)
aic.dt <- aic.dt[abs(aic.dt$Estimate) > .1, ]
rownames(aic.dt) <- gsub("\\(Intercept\\)", "Intercept", rownames(aic.dt))
aic.dt$significance <- -log10(aic.dt$"Pr(>|t|)")
aic.dt$TF <- rep("NO", nrow(aic.dt))
aic.dt$TF[grep("ctcf", rownames(aic.dt))] <- "YES"
aic.dt$TF[grep("pol2", rownames(aic.dt))] <- "YES"
aic.dt$TF <- factor(aic.dt$TF, levels = c("YES", "NO"))
pdf("~/Dropbox/Network/R1.fig/effectAIC.pdf")
ggplot(aic.dt, aes(x = rownames(aic.dt), y = Estimate)) +  
  geom_point(aes(color = TF, size = significance)) + 
  # scale_y_continuous(lim = c(-.2, .2)) +
  theme_bw() + xlab("") + ylab("Effect") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold")) + 
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"), legend.key = element_blank()) 
dev.off()
