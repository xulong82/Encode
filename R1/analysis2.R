library(parallel)

rm(list = ls())
load("~/Dropbox/Network/R1/esb4Df.rdt")  # ESB4
load("~/Dropbox/Network/R1/esb4Analysis.rdt")  # ESB4

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

geno <- log2(as.matrix(data1[, -c(1:11)]) + 1)  # log2
geno.sd <- apply(geno, 2, sd)
geno <- geno[, geno.sd > 0]

data1 <- as.data.frame(pair)
data1$label <- apply(pairs, 2, function (x) {paste(x, collapse = " X ")})
data1$q.value <- p.adjust(data1$p.value, method = "fdr")
data1 <- data1[data1$q.value < .01, ]
data1 <- data1[data1$logLik > quantile(data1$logLik, .95), ]
data1 <- data1[abs(data1$interaction) > 3, ]
data1$mark1 <- gsub("\\..*", "", data1$label)
data1$mark2 <- gsub("^.*X ", "", data1$label)
data1$mark2 <- gsub("\\..*", "", data1$mark2)
x = t(apply(data1[, 2:3], 1, sort))
data1$pair <- paste(x[, 1], x[, 2], sep = " X ")
data1 <- data1[data1$mark1 != data1$mark2, ]

#--- INTERACTION PLOT
#--- select the features
features <- c("k4me1.tx.intensity", "k4me3.tx.intensity")
features <- c("k9ac.tx.intensity", "pol2.tx.intensity")

data3 <- data.frame(rna = rna2, geno[, features])
summary(data3[, 1])
summary(data3[, 2])
summary(data3[, 3])

pdf("~/Dropbox/Network/R1.fig/pairs.pdf")
pairs(data3)
dev.off()

data3$f1 <- data3$f2 <- rep("q2", nrow(data3))

data3$f1[data3[, 2] == 0] <- "q1"
data3$f1[data3[, 2] > 0] <- "q3"
data3$f2[data3[, 3] == 0] <- "q1"
data3$f2[data3[, 3] > 0] <- "q3"

data3$f1[data3[, 2] < mean(data3[, 2])] <- "q1"
data3$f1[data3[, 2] > mean(data3[, 2])] <- "q3"
data3$f2[data3[, 3] < mean(data3[, 3])] <- "q1"
data3$f2[data3[, 3] > mean(data3[, 3])] <- "q3"

data3$group <- paste(data3$f1, data3$f2, sep = "/")
data4 <- data3

data3$f1[data3[, 2] < quantile(data3[, 2], .25)] <- "q1"
data3$f1[data3[, 2] > quantile(data3[, 2], .75)] <- "q3"
data3$f2[data3[, 3] < quantile(data3[, 3], .25)] <- "q1"
data3$f2[data3[, 3] > quantile(data3[, 3], .75)] <- "q3"

data3$group <- paste(data3$f1, data3$f2, sep = "/")
data4 <- data3[-grep("q2", data3$group), ]
# data4$group <- gsub("q1", paste(features[1], "q1", sep = "-"), data4$group)
# data4$group <- gsub("q3", paste(features[2], "q3", sep = "-"), data4$group)
table(data4$group)

pdf("~/Dropbox/Network/R1.fig/q1q3plot1.pdf", width = 6)
ggplot(data4, aes(x = group, y = rna, fill = group)) + geom_boxplot() +
# geom_hline(yintercept = 65, color = "darkred", size = 2, linetype = 2) + 
  theme_bw() + xlab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

model1 <- lm(as.formula(paste("rna ~ ", paste(colnames(data3)[2], collapse= "+"))), data3)  # additive
model2 <- lm(as.formula(paste("rna ~ ", paste(colnames(data3)[3], collapse= "+"))), data3)  # additive
model3 <- lm(as.formula(paste("rna ~ ", paste(colnames(data3)[2:3], collapse= "+"))), data3)  # additive
model4 <- lm(as.formula(paste("rna ~ ", paste(colnames(data3)[2:3], collapse= "*"))), data3)  # complete

model1.dt <- data.frame(feature = names(model1$coefficients), effect = model1$coefficients)
model1.dt$model <- rep("model1", nrow(model1.dt))
model2.dt <- data.frame(feature = names(model2$coefficients), effect = model2$coefficients)
model2.dt$model <- rep("model2", nrow(model2.dt))
model3.dt <- data.frame(feature = names(model3$coefficients), effect = model3$coefficients)
model3.dt$model <- rep("model3", nrow(model3.dt))
model4.dt <- data.frame(feature = names(model4$coefficients), effect = model4$coefficients)
model4.dt$model <- rep("model4", nrow(model4.dt))
model.dt <- rbind(model1.dt, model2.dt, model3.dt, model4.dt)

model.dt$feature <- gsub("\\(Intercept\\)", "Intercept", model.dt$feature)
model.dt$position <- as.numeric(gsub("model", "", model.dt$model))
model.dt$feature <- factor(model.dt$feature, levels = unique(model.dt$feature))

pdf("~/Dropbox/Network/R1.fig/model1.pdf", width = 6)
ggplot(model.dt, aes(x = position, y = effect, color = feature) ) +  
  geom_point(size = 5) + geom_line(linetype = 3, width = 2) +
  theme_bw() + xlab("") + ylab("Effect") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()
