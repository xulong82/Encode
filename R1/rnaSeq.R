library(amap)
library(ape)
library(ggplot2)

rm(list = ls())
#--- Encode: RNA-seq
name1 <- list.files(path = "~/Dropbox/Network/RSEM2", pattern = "*.genes.results")
name2 <- gsub(".genes.results", "", name1)
name3 <- unique(gsub("Rep.*", "", name1))  # sample name
for (i in 1:length(name1)) {
  filepath <- file.path("~/Dropbox/Network/RSEM2", name1[i])
  cat(i, "/", length(name1), name1[i], "\n")
  assign(name2[i], read.delim(filepath, stringsAsFactors = F)[, -2])
}

rna <- matrix(nrow = nrow(get(name2[1])), ncol = length(name2),  # rna: TPM 
              dimnames = list(get(name2[1])$gene_id, name2))
for (i in 1:length(name2)) rna[, i] <- get(name2[i])$TPM

ens2sym <- read.delim("~/Dropbox/X/ensembl2symbol.map", header = F, row.names = 1, stringsAsFactors = F)
rna <- as.data.frame(rna)
rna$symbol <- ens2sym[rownames(rna), ]
rna <- aggregate(. ~ symbol, rna, sum)
rna <- data.frame(row.names = rna$symbol, rna[, 2:ncol(rna)])

rna1 <- rna[rowSums(rna) > ncol(rna), ]  # miminum
rna1 <- rna1[apply(rna1, 1, function (x) {max(x) > 1e2}), ]  # maximal
rna1 <- log2(rna1 + 1)  # log2

tissue <- gsub("Rep[12]", "", colnames(rna1))
aov.pval <- rep(NA, nrow(rna1))
for (i in 1:nrow(rna1)) {  # ANOVA
  if (i %% 2e3 == 0) cat(i, "/", nrow(rna1), "\n")
  dat1 <- data.frame(tpm = as.matrix(rna1)[i, ], tissue)
  aov1 = aov(tpm ~ tissue, data = dat1)
  aov.pval[i] <- min(summary(aov1)[[1]][["Pr(>F)"]], na.rm = T)
}
rna2 <- rna1[aov.pval < 0.01, ]  # anova
rna2 <- rna2[apply(rna2, 1, function (x) {max(x) - min(x) > 2}), ]  # 4 fold changes

tissue1 <- unique(tissue)
rna2.mean <- matrix(nrow = nrow(rna2), ncol = length(tissue1), dimnames = list(rownames(rna2), tissue1))
for (i in 1:length(tissue1)) rna2.mean[, i] <- rowMeans(rna2[, grep(tissue1[i], colnames(rna2))])
rna.df <- rna2.mean
save(rna.df, file = "~/Dropbox/Network/R1/rna.rdt")

rna3 <- rna2[rna2.mean[, "Testis"] > apply(rna2.mean, 1, function(x) {quantile(x, 0.75)}), ]  # 0.75 quantile
rna3.mean <- rna2.mean[rownames(rna3), ]

genes <- rownames(rna3)
write(genes, file = "~/Dropbox/Network/Genes/testisGenes.txt")
save.image("~/Dropbox/Network/R1/genes.rdt")

rsem1 <- read.delim("~/Dropbox/Network/RSEM2/rsem1.txt", header = F, stringsAsFactors = F)
rsem2 <- read.delim("~/Dropbox/Network/RSEM2/rsem2.txt", header = F, stringsAsFactors = F)
rsem1 <- gsub("^.*: ", "", rsem1$V1)
rsem2 <- gsub("FASTQ2\\/", "", rsem2$V1)
rsem2 <- gsub("_cut\\.fastq", "", rsem2)
rsem.count <- as.numeric(gsub(" .*", "", rsem1))
rsem.percent <- gsub("^.*\\(", "", rsem1)
rsem.percent <- as.numeric(gsub("\\%.*", "", rsem.percent))

rsem <- data.frame(sample = rsem2, rsem.count, rsem.percent)
rsem$rsem.count <- rsem$rsem.count * 1e-6
rsem$type <- gsub("Rep[12]", "", rsem$sample)
# pdf("~/Dropbox/Network/R1.fig/rsem1.pdf", width = 6)
pdf("~/Dropbox/Network/R1.fig/rsem2.pdf", width = 6)
# ggplot(rsem, aes(x = sample, y = rsem.count)) +  
ggplot(rsem, aes(x = sample, y = rsem.percent)) +  
  geom_point(aes(color = type), size = 3) + 
  theme_bw() + xlab("") + ylab("Aligned read count (M)") +
  theme_bw() + xlab("") + ylab("Alignment rate (%)") +
  # theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, color = "gray30", face = "bold"),
        axis.text.y = element_text(size = 10, color = "gray30", face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

# heatmap(cor(rna, method = "pearson"))
hc1 <- hcluster(t(rna1), method = "pearson", link = "average")  # clust samples
hc2 <- hcluster(rna1, method = "pearson", link = "average")  # clust genes
module <- cutree(hc2, h = .25)  # pho > .75
table(names(module) == rownames(rna1)) 
pdf("~/Dropbox/Network/R1.fig/phylo1.pdf")
# plot(as.phylo(hc1))
plot(as.phylo(hc1), type = "unrooted", lab4ut = "axial")
dev.off()
