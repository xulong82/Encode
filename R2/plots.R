#---
mark <- h3k4me1

scan1.h3k4me1$mark <- rep("h3k4me1", nrow(scan1.h3k4me1))
scan1.h3k4me3$mark <- rep("h3k4me3", nrow(scan1.h3k4me3))
scan1.h3k27ac$mark <- rep("h3k27ac", nrow(scan1.h3k27ac))
scan1.h3k27me3$mark <- rep("h3k27me3", nrow(scan1.h3k27me3))

scan1.h3k4me1$pos 

id.gene1 <- ids.rna.mod[[1]]

scan1.h3k4me1.2$gene <- as.character(scan1.h3k4me1.2$gene)
scan1.h3k4me1.2$epi <- as.character(scan1.h3k4me1.2$epi)
scan1.h3k4me1.2$gene.pos <- ens.ucsc[scan1.h3k4me1.2$gene, "pos.Mbp"]
scan1.h3k4me1.2$epi.pos <- rowMeans(epi.mks[scan1.h3k4me1.2$epi, ]) * 2e2 * 1e-6

scan1.mmark.2$gene <- as.character(scan1.mmark.2$gene)
scan1.mmark.2$epi <- as.character(scan1.mmark.2$epi)
scan1.mmark.2$gene.pos <- ens.ucsc[scan1.mmark.2$gene, "pos.Mbp"]
scan1.mmark.2$epi.pos <- rowMeans(epi.mks[scan1.mmark.2$epi, ]) * 2e2 * 1e-6

pdf(paste(paste("~/Dropbox/Network/R2.fig/map", mark, sep = "."), "pdf", sep = "."))
pdf("~/Dropbox/Network/R2.fig/map.mmark.pdf")
ggplot(scan1.mmark.2, aes(x = gene.pos, y = epi.pos)) + 
  geom_point(color = "blue", size = .2) + 
# geom_point(aes(color = logLik), size = .2) + 
# geom_point(color = "red", aes(size = logLik)) + 
  theme_bw() + 
  xlab("Gene") + ylab("Epispots") +
  scale_x_continuous(breaks = chrmid.Mbp, labels = names(chrlen2)) +
  scale_y_continuous(breaks = chrmid.Mbp, labels = names(chrlen2)) +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()
#---
idx1.epi$chrom <- pos200tochrom[idx1.epi$start, "chrom"]
idx1.epi$start1 <- idx1.epi$start - c(0, chrlen3)[idx1.epi$chrom]
idx1.epi$end1 <- idx1.epi$end- c(0, chrlen3)[idx1.epi$chrom]
idx1.epi$center <- (idx1.epi$start1 + idx1.epi$end1) /2

pdf("~/Dropbox/Network/R2.fig/dis.spot1.pdf", width = 12)
ggplot(idx1.epi, aes(x = center, y = 1)) +
  geom_point(color = "black", size = .1) +
  geom_vline(aes(xintercept = center)) +
  facet_grid(chrom ~ .) +
  theme_bw() +
  xlab("Genome position (200 bp unit)") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_blank())
dev.off()
