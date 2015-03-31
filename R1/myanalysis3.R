# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Encode data integration - network decomposition
# Rev: March 13, 2014

source("~/Dropbox/Encode/R/myfunctions.R")

refgene.tss.b <- myBinary(refgene.tss)

geno.k4me1 <- myBinary(testis.h3k4me1)[refgene.tss.b == 1]  # A
geno.k4me3 <- myBinary(testis.h3k4me3)[refgene.tss.b == 1]  # B
geno.k27me3 <- myBinary(testis.h3k27me3)[refgene.tss.b == 1]  # C

geno.WT <- which(geno.k4me1 == 0 & geno.k4me3 == 0 & geno.k27me3 == 0)
geno.A <- which(geno.k4me1 == 1 & geno.k4me3 == 0 & geno.k27me3 == 0)
geno.B <- which(geno.k4me1 == 0 & geno.k4me3 == 1 & geno.k27me3 == 0)
geno.C <- which(geno.k4me1 == 0 & geno.k4me3 == 0 & geno.k27me3 == 1)
geno.AB <- which(geno.k4me1 == 1 & geno.k4me3 == 1 & geno.k27me3 == 0)
geno.AC <- which(geno.k4me1 == 1 & geno.k4me3 == 0 & geno.k27me3 == 1)
geno.BC <- which(geno.k4me1 == 0 & geno.k4me3 == 1 & geno.k27me3 == 1)

pheno.rna <- myBinary(testis.rna)[refgene.tss.b == 1]

table(pheno.rna[geno.AB])
table(pheno.rna[geno.AC])
table(pheno.rna[geno.BC])
