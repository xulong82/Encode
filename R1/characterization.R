rm(list = ls())
load("~/Dropbox/Network/R1/encodex.rdt")

load("~/Dropbox/Network/R1/esb4Df.rdt")  # ESB4
load("~/Dropbox/Network/R1/melDf.rdt")  # MEL
load("~/Dropbox/Network/R1/liverDf.rdt")  # Liver

marks.bi <- marks.bi[, -10]  # delete p300

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

library(parallel)
y <- mclapply(mylist, function (x) {apply(x, 2, Peak)}, mc.cores = 5)
y1 <- matrix(unlist(y), nrow = 5, byrow = T)

number <- rbind(PM2k = y[[3]], Exon = y[[4]], Intron = y[[5]], Intergene = y[[1]] - y[[2]])
number <- number * 1e-4  #
  
# save(y, length, number, file = "~/Dropbox/Network/R1/characterizationEsb4.rdt")
# save(y, length, number, file = "~/Dropbox/Network/R1/characterizationMel.rdt")
save(y, length, number, file = "~/Dropbox/Network/R1/characterizationLiver.rdt")

#--- graphs
load("~/Dropbox/Network/R1/characterizationEsb4.rdt")
load("~/Dropbox/Network/R1/characterizationMel.rdt")
load("~/Dropbox/Network/R1/characterizationLiver.rdt")
library(ggplot2)
#----------------------------------------------------------------------------------------------
bar.dt1 <- data.frame(value = c(as.vector(length)), 
                      mark = rep(colnames(length), each = nrow(length)),
                      feature = rep(rownames(length), ncol(length)))
bar.dt1$order <- reorder(bar.dt1$mark, bar.dt1$value)
bar.dt2 <- data.frame(value = c(as.vector(number)), 
                      mark = rep(colnames(number), each = nrow(number)),
                      feature = rep(rownames(number), ncol(number)))
bar.dt2$order <- reorder(bar.dt2$mark, bar.dt2$value)

pdf("~/Dropbox/Network/R1.fig/char1.pdf", width = 4)
# pdf("~/Dropbox/Network/R1.fig/char2.pdf", width = 4)
ggplot(bar.dt1, aes(x = order, y = value, fill = feature)) +  # bar.dt1, bar.dt2
# ggplot(bar.dt2, aes(x = order, y = value, fill = feature)) +  # bar.dt1, bar.dt2
  geom_bar(stat = "identity") +
# geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  xlab("") + ylab("Genome %") +  # bar.dt1
# xlab("") + ylab("Peak number x 1e-4") +  # bar.dt2
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

bar.dt3$order <- reorder(bar.dt3$mark, bar.dt3$value)

pdf("~/Dropbox/Network/R1.fig/char3.pdf", width = 10)
ggplot(bar.dt3, aes(x = order, y = value, fill = type1)) +  # bar.dt1, bar.dt2
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

#--- peaks over the genome
setwd("~/Dropbox/Network/ChipSeq/mm10")
h3k36me3 <- Colname(read.delim("HistoneEsb4H3k36me3ME0.bed", header =F, stringsAsFactors = F))
h3k36me3$pos <- (h3k36me3$chromStart + h3k36me3$chromEnd) * 1e-6 / 2
h3k36me3$length <- (h3k36me3$chromEnd - h3k36me3$chromStart) * 1e-6
h3k36me3$chrom <- factor(h3k36me3$chrom, levels = c(paste("chr", c(1:19, "X", "Y"), sep = "")))
pdf(file = "~/Dropbox/Network/R1.fig/broadpeak1.pdf", width = 15, height = 10)
ggplot(h3k36me3, aes(x = pos, y = 1)) +
  geom_point(aes(color = length, size = length)) +
  # geom_point(color = "red", size = .1) +
  # geom_vline(aes(xintercept = pos)) +
  facet_grid(chrom ~ .) +
  theme_bw() +
  xlab("Genome (Mbp)") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, face = "bold", color = "grey30"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank()) +
  theme(legend.position = "none")
dev.off()
