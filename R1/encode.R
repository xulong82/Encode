# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - import the encode data
# Rev: March 13, 2014

#--- ChIP-seq broadPeak data from Encode
setwd("~/Dropbox/Encode/ChIP-seq/broadPeak")

brain.h3k4me1 <- read.delim("wgEncodeLicrHistoneWbrainH3k04me1UE14halfC57bl6StdPk.broadPeak", header = F)
brain.h3k4me3 <- read.delim("wgEncodeLicrHistoneWbrainH3k04me3UE14halfC57bl6StdPk.broadPeak", header = F)
brain.h3k27me3 <- read.delim("wgEncodeLicrHistoneWbrainH3k27me3ME14halfC57bl6StdPk.broadPeak", header = F)
brain.h3k27ac <- read.delim("wgEncodeLicrHistoneWbrainH3k27acUE14halfC57bl6StdPk.broadPeak", header = F)
colnames(brain.h3k4me1)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(brain.h3k4me3)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(brain.h3k27me3)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(brain.h3k27ac)[1:3] <- c("chrom", "chromStart", "chromEnd")

spleen.h3k4me1 <- read.delim("wgEncodeLicrHistoneSpleenH3k4me1MAdult8wksC57bl6StdPk.broadPeak", header = F)
spleen.h3k4me3 <- read.delim("wgEncodeLicrHistoneSpleenH3k4me3MAdult8wksC57bl6StdPk.broadPeak", header = F)
spleen.h3k27me3 <- read.delim("wgEncodeLicrHistoneSpleenH3k27me3MAdlt8wC57bl6StdPk.broadPeak", header = F)
spleen.h3k27ac <- read.delim("wgEncodeLicrHistoneSpleenH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F) 
colnames(spleen.h3k4me1)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(spleen.h3k4me3)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(spleen.h3k27me3)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(spleen.h3k27ac)[1:3] <- c("chrom", "chromStart", "chromEnd")

testis.h3k4me1 <- read.delim("wgEncodeLicrHistoneTestisH3k04me1MAdult8wksC57bl6StdPk.broadPeak", header = F)
testis.h3k4me3 <- read.delim("wgEncodeLicrHistoneTestisH3k04me3MAdult8wksC57bl6StdPk.broadPeak", header = F)
testis.h3k27me3 <- read.delim("wgEncodeLicrHistoneTestisH3k27me3MAdlt8wC57bl6StdPk.broadPeak", header = F)
testis.h3k27ac <- read.delim("wgEncodeLicrHistoneTestisH3k27acMAdult8wksC57bl6StdPk.broadPeak", header = F)
colnames(testis.h3k4me1)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(testis.h3k4me3)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(testis.h3k27me3)[1:3] <- c("chrom", "chromStart", "chromEnd")
colnames(testis.h3k27ac)[1:3] <- c("chrom", "chromStart", "chromEnd")

save.image("~/Dropbox/Encode/R/chipseq.rdt")
load("~/Dropbox/Encode/R/chipseq.rdt")

#--- RNA-seq gene expression data from Encode
setwd("~/Dropbox/Encode/RNA-seq")
brain.rna <- read.delim("wgEncodeCshlLongRnaSeqWbrainE14halfGeneEnsV65IAcuff-Rev.gtf", header = F)
spleen.rna <- read.delim("wgEncodeCshlLongRnaSeqSpleenAdult8wksGeneEnsV65IAcuff-Rev.gtf", header = F)
testis.rna <- read.delim("wgEncodeCshlLongRnaSeqTestisAdult8wksGeneEnsV65IAcuff-Rev.gtf", header = F)

brain.rna <- data.frame(chrom      = brain.rna$V1, 
                        chromStart = brain.rna$V4, 
                        chromEnd   = brain.rna$V5, 
                        gene       = brain.rna$V9, 
                        RPKM1      = brain.rna$V11,
                        RPKM2      = brain.rna$V12,
                        RPKM       = rowMeans(brain.rna[11:12]))
spleen.rna <- data.frame(chrom      = spleen.rna$V1, 
                         chromStart = spleen.rna$V4, 
                         chromEnd   = spleen.rna$V5, 
                         gene       = spleen.rna$V9, 
                         RPKM1      = spleen.rna$V11,
                         RPKM2      = spleen.rna$V12,
                         RPKM       = rowMeans(spleen.rna[11:12]))
testis.rna <- data.frame(chrom      = testis.rna$V1, 
                         chromStart = testis.rna$V4, 
                         chromEnd   = testis.rna$V5, 
                         gene       = testis.rna$V9, 
                         RPKM1      = testis.rna$V11,
                         RPKM2      = testis.rna$V12,
                         RPKM       = rowMeans(testis.rna[11:12]))
brain.rna <- unique(brain.rna)

ensembl2symbol.map <- read.delim("~/Dropbox/Encode/ensembl2symbol.map", sep = "\t", header = F)
brain.rna <- merge(brain.rna, ensembl2symbol.map, by.x = "gene", by.y = "V1", all.x = T)
spleen.rna <- merge(spleen.rna, ensembl2symbol.map, by.x = "gene", by.y = "V1", all.x = T)
testis.rna <- merge(testis.rna, ensembl2symbol.map, by.x = "gene", by.y = "V1", all.x = T)

brain.rna <- brain.rna[brain.rna$RPKM > 2, ]
spleen.rna <- spleen.rna[spleen.rna$RPKM > 2, ]
testis.rna <- testis.rna[testis.rna$RPKM > 2, ]

brain.rna <- brain.rna[brain.rna$chrom != "chrMT", ]
spleen.rna <- spleen.rna[spleen.rna$chrom != "chrMT", ]
testis.rna <- testis.rna[testis.rna$chrom != "chrMT", ]

brain.rna <- data.frame(chrom      = brain.rna$chrom,  # Embryonic 14.5 days 
                           chromStart = brain.rna$chromStart,
                           chromEnd   = brain.rna$chromEnd,
                           gene       = brain.rna$V2,
                           RPKM       = brain.rna$RPKM)
spleen.rna <- data.frame(chrom      = spleen.rna$chrom,  # Adult 8 weeks 
                            chromStart = spleen.rna$chromStart,
                            chromEnd   = spleen.rna$chromEnd,
                            gene       = spleen.rna$V2,
                            RPKM       = spleen.rna$RPKM)
testis.rna <- data.frame(chrom      = testis.rna$chrom,  # Adult 8 weeks 
                            chromStart = testis.rna$chromStart,
                            chromEnd   = testis.rna$chromEnd,
                            gene       = testis.rna$V2,
                            RPKM       = testis.rna$RPKM)

brain.rna <- aggregate(RPKM ~ ., data = brain.rna, max)
spleen.rna <- aggregate(RPKM ~ ., data = spleen.rna, max)
testis.rna <- aggregate(RPKM ~ ., data = testis.rna, max)

save(brain.rna, spleen.rna, testis.rna, file = "~/Dropbox/Encode/R/rna.RData")  # RPKM > 2
save.image("~/Dropbox/Encode/R/encode.20140313.RData")
