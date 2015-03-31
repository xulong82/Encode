# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network
# Rev: August 10, 2014

chr <- read.delim("~/Dropbox/X/chromInfo.txt", header = F)  # chromosomes length
chr.200 <- as.numeric(chr$V2) %/% 200 + 1

chrlen1 <- cumsum(as.numeric(chr$V2))  # bp unit
chrlen2 <- chrlen1 * 1e-6  # Mbp unit
chrlen3 <- chrlen1 %/% 200 + 1  #  200 bp unit 
names(chrlen1) <- 1:21
names(chrlen2) <- 1:21
names(chrlen3) <- 1:21
genome.200 <- chrlen3[length(chrlen3)]

chrmid.bp <- diff(c(0, chrlen1)) * 0.5 + c(0, chrlen1[-length(chrlen1)])  # Mbp unit
chrmid.Mbp <- diff(c(0, chrlen2)) * 0.5 + c(0, chrlen2[-length(chrlen2)])  # Mbp unit

ens.ucsc <- read.delim("~/Dropbox/X/myRefGene.tx", stringsAsFactors = F)
ens.ucsc <- ens.ucsc[!duplicated(ens.ucsc$name2), ]
ens.ucsc <- data.frame(row.names  = ens.ucsc$name2, 
                       chrom      = ens.ucsc$chrom,
                       chromStart = ens.ucsc$txStart,
                       chromEnd   = ens.ucsc$txEnd)
ens.ucsc$chrom <- as.character(ens.ucsc$chrom)
ens.ucsc$chrom <- gsub("chrX", "chr20", ens.ucsc$chrom)
ens.ucsc$chrom <- gsub("chrY", "chr21", ens.ucsc$chrom)
ens.ucsc$chrom <- gsub("chr", "", ens.ucsc$chrom)
ens.ucsc$chrom <- as.numeric(ens.ucsc$chrom)

ens.ucsc$pos.Mbp <- c(0, chrlen2)[ens.ucsc$chrom] + rowMeans(ens.ucsc[, 2:3]) * 1e-6

pos200tochrom <- data.frame(pos.200 = 1:genome.200, chrom = rep(0, genome.200))
for (i in 1:length(chrlen1)) {
  start <- c(0, chrlen3)[i] + 1
  end <- chrlen3[i]
# cat("chromosome", i, start, end, "\n")
  pos200tochrom$chrom[start:end] = i
}
