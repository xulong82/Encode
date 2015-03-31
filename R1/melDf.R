# Copyright: Xulong Wang (xulong.wang@jax.org)
# Rev: August 28, 2014

rm(list = ls())
# source("~/Dropbox/Network/R1/encodex.R")
load("~/Dropbox/Network/R1/encodex.rdt")
load("~/Dropbox/Network/R1/rnaSeq.rdt")

#--- ENCODE broadPeak in mm9 lifted over to BED in mm10 assembly
setwd("~/Dropbox/Network/ChipSeq/mm10")
h3k4me1 <- Colname(read.delim("HistoneMelH3k04me1MImmortal.bed", header = F, stringsAsFactors = F))
h3k4me3 <- Colname(read.delim("HistoneMelH3k04me3MImmortal.bed", header = F, stringsAsFactors = F))
h3k9ac <- Colname(read.delim("HistoneMelH3k09acMImmortal.bed", header = F, stringsAsFactors = F))
# h3k9me3 <- Colname(read.delim("", header = F, stringsAsFactors = F))
h3k27ac <- Colname(read.delim("HistoneMelH3k27acMImmortal.bed", header = F, stringsAsFactors = F))
h3k27me3 <- Colname(read.delim("HistoneMelH3k27me3MImmortal.bed", header = F, stringsAsFactors = F))
h3k36me3 <- Colname(read.delim("HistoneMelH3k36me3MImmortal.bed", header = F, stringsAsFactors = F))
h3k79me2 <- Colname(read.delim("HistoneMelH3k79me2MImmortal.bed", header = F, stringsAsFactors = F))
ctcf <- Colname(read.delim("TfbsMelCtcfMImmortal.bed", header = F, stringsAsFactors = F))
pol2 <- Colname(read.delim("TfbsMelPol2MImmortal.bed", header = F, stringsAsFactors = F))
# p300 <- Colname(read.delim("", header = F, stringsAsFactors = F))

#--- Segment ENCODE broadPeak
marks.bi <- cbind(Binary(h3k4me1), Binary(h3k4me3), Binary(h3k9ac), # Binary(h3k9me3),
  Binary(h3k27ac), Binary(h3k27me3), Binary(h3k36me3), Binary(h3k79me2), Binary(ctcf), Binary(pol2))
marks.si <- cbind(Signal(h3k4me1), Signal(h3k4me3), Signal(h3k9ac), # Signal(h3k9me3),
  Signal(h3k27ac), Signal(h3k27me3), Signal(h3k36me3), Signal(h3k79me2), Signal(ctcf), Signal(pol2))
colnames(marks.si) <- colnames(marks.bi) <- c("k4me1", "k4me3", "k9ac", "k27ac", "k27me3", "k36me3", "k79me2", "ctcf", "pol2")
rownames(marks.si) <- rownames(marks.bi) <- 1:nrow(marks.bi)

#--- define enhancer pseduo-mark by mark combinations: active/inactive
# active.en <- as.numeric(as.logical(marks.bi[, "k4me1"]) & 
#   as.logical(marks.bi[, "k27ac"]) & !as.logical(marks.bi[, "k27me3"]))
# inactive.en <- as.numeric(as.logical(marks.bi[, "k4me1"]) & 
#   !as.logical(marks.bi[, "k27ac"]) & as.logical(marks.bi[, "k27me3"]))

#--- epigenetic features in genetic regions
refgene.marks.si <- refgene.marks.bi <- array(dim = c(ncol(refgene.bi), nrow(marks.bi), ncol(marks.bi)),
                                              dimnames = list(colnames(refgene.bi), rownames(marks.bi), colnames(marks.bi)))
for (i in 1:ncol(refgene.bi)) {
  cat(i, "in", ncol(refgene.bi), "\n")
  refgene.marks.bi[i, , ] <- refgene.bi[, i] * marks.bi
  refgene.marks.si[i, , ] <- refgene.bi[, i] * marks.si
}

#--- refine the refgene annotation entries
rna.df <- rna.df[(rownames(rna.df) %in% rownames(refgene.tx)), ]
refgene.tx <- refgene.tx[rownames(rna.df), ]
refgene.tss1k <- refgene.tss1k[rownames(rna.df), ]
refgene.pm2k <- refgene.pm2k[rownames(rna.df), ]
refgene.pm5k <- refgene.pm5k[rownames(rna.df), ]
refgene.en <- refgene.en[rownames(rna.df), ]

#--- lengths of genetic features
tx.len <- (refgene.tx$chromEnd - refgene.tx$chromStart) %/% chrom.unit + 1
ex.len <- Geno(refgene.tx, refgene.bi[, "ex"], "length")
it.len <- Geno(refgene.tx, refgene.bi[, "it"], "length")
u5.len <- Geno(refgene.tx, refgene.bi[, "u5"], "length")
u3.len <- Geno(refgene.tx, refgene.bi[, "u3"], "length")

#--- choose genes by transcript length
idx <- unique(c(which(tx.len < 2), which(ex.len < 2), which(it.len < 2), which(u5.len == 0), which(u3.len == 0)))
refgene.tx <- refgene.tx[-idx, ]
refgene.tss1k <- refgene.tss1k[-idx, ]
refgene.pm2k <- refgene.pm2k[-idx, ]
refgene.pm5k <- refgene.pm5k[-idx, ]
refgene.en <- refgene.en[-idx, ]

cat(date(), "--- I: total length in given feature \n")
tx.length <- apply(refgene.marks.bi["tx", , ], 2, function (x) {Geno(refgene.tx, x, "length")})
tss1k.length <- apply(refgene.marks.bi["tss1k", , ], 2, function (x) {Geno(refgene.tss1k, x, "length")})
pm2k.length <- apply(refgene.marks.bi["pm2k", , ], 2, function (x) {Geno(refgene.pm2k, x, "length")})
pm5k.length <- apply(refgene.marks.bi["pm5k", , ], 2, function (x) {Geno(refgene.pm5k, x, "length")})
ex.length <- apply(refgene.marks.bi["ex", , ], 2, function (x) {Geno(refgene.tx, x, "length")})
it.length <- apply(refgene.marks.bi["it", , ], 2, function (x) {Geno(refgene.tx, x, "length")})
en.length <- apply(refgene.marks.bi["en", , ], 2, function (x) {Geno(refgene.en, x, "length")})
u5.length <- apply(refgene.marks.bi["u5", , ], 2, function (x) {Geno(refgene.tx, x, "length")})
u3.length <- apply(refgene.marks.bi["u3", , ], 2, function (x) {Geno(refgene.tx, x, "length")})

colnames(tx.length) <- paste(colnames(tx.length), "tx.length", sep = ".")
colnames(tss1k.length) <- paste(colnames(tss1k.length), "tss1k.length", sep = ".")
colnames(pm2k.length) <- paste(colnames(pm2k.length), "pm2k.length", sep = ".")
colnames(pm5k.length) <- paste(colnames(pm5k.length), "pm5k.length", sep = ".")
colnames(ex.length) <- paste(colnames(ex.length), "ex.length", sep = ".")
colnames(it.length) <- paste(colnames(it.length), "it.length", sep = ".")
colnames(en.length) <- paste(colnames(en.length), "en.length", sep = ".")
colnames(u5.length) <- paste(colnames(u5.length), "u5.length", sep = ".")
colnames(u3.length) <- paste(colnames(u3.length), "u3.length", sep = ".")

cat(date(), "--- II: signal peak in given feature \n")
tx.signal <- apply(refgene.marks.si["tx", , ], 2, function (x) {Geno(refgene.tx, x, "signal")})
tss1k.signal <- apply(refgene.marks.si["tss1k", , ], 2, function (x) {Geno(refgene.tss1k, x, "signal")})
pm2k.signal <- apply(refgene.marks.si["pm2k", , ], 2, function (x) {Geno(refgene.pm2k, x, "signal")})
pm5k.signal <- apply(refgene.marks.si["pm5k", , ], 2, function (x) {Geno(refgene.pm5k, x, "signal")})
ex.signal <- apply(refgene.marks.si["ex", , ], 2, function (x) {Geno(refgene.tx, x, "signal")})
it.signal <- apply(refgene.marks.si["it", , ], 2, function (x) {Geno(refgene.tx, x, "signal")})
en.signal <- apply(refgene.marks.si["en", , ], 2, function (x) {Geno(refgene.en, x, "signal")})
u5.signal <- apply(refgene.marks.si["u5", , ], 2, function (x) {Geno(refgene.tx, x, "signal")})
u3.signal <- apply(refgene.marks.si["u3", , ], 2, function (x) {Geno(refgene.tx, x, "signal")})

colnames(tx.signal) <- paste(colnames(tx.signal), "tx.signal", sep = ".")
colnames(tss1k.signal) <- paste(colnames(tss1k.signal), "tss1k.signal", sep = ".")
colnames(pm2k.signal) <- paste(colnames(pm2k.signal), "pm2k.signal", sep = ".")
colnames(pm5k.signal) <- paste(colnames(pm5k.signal), "pm5k.signal", sep = ".")
colnames(ex.signal) <- paste(colnames(ex.signal), "ex.signal", sep = ".")
colnames(it.signal) <- paste(colnames(it.signal), "it.signal", sep = ".")
colnames(en.signal) <- paste(colnames(en.signal), "en.signal", sep = ".")
colnames(u5.signal) <- paste(colnames(u5.signal), "u5.signal", sep = ".")
colnames(u3.signal) <- paste(colnames(u3.signal), "u3.signal", sep = ".")

cat(date(), "--- III: mark intensity in given feature \n")
tx.intensity <- tx.length / tx.len[-idx]
ex.intensity <- ex.length / ex.len[-idx]
it.intensity <- it.length / it.len[-idx]
u5.intensity <- u5.length / u5.len[-idx]
u3.intensity <- u3.length / u3.len[-idx]
colnames(tx.intensity) <- gsub("length", "intensity", colnames(tx.intensity))
colnames(ex.intensity) <- gsub("length", "intensity", colnames(ex.intensity))
colnames(it.intensity) <- gsub("length", "intensity", colnames(it.intensity))
colnames(u5.intensity) <- gsub("length", "intensity", colnames(u5.intensity))
colnames(u3.intensity) <- gsub("length", "intensity", colnames(u3.intensity))

cat(date(), "--- IV: number of peaks in given feature \n")
tx.peak <- apply(refgene.marks.bi["tx", , ], 2, function (x) {Geno(refgene.tx, x, "peak")})
tss1k.peak <- apply(refgene.marks.bi["tss1k", , ], 2, function (x) {Geno(refgene.tss1k, x, "peak")})
pm2k.peak <- apply(refgene.marks.bi["pm2k", , ], 2, function (x) {Geno(refgene.pm2k, x, "peak")})
pm5k.peak <- apply(refgene.marks.bi["pm5k", , ], 2, function (x) {Geno(refgene.pm5k, x, "peak")})
ex.peak <- apply(refgene.marks.bi["ex", , ], 2, function (x) {Geno(refgene.tx, x, "peak")})
it.peak <- apply(refgene.marks.bi["it", , ], 2, function (x) {Geno(refgene.tx, x, "peak")})
en.peak <- apply(refgene.marks.bi["en", , ], 2, function (x) {Geno(refgene.en, x, "peak")})
u5.peak <- apply(refgene.marks.bi["u5", , ], 2, function (x) {Geno(refgene.tx, x, "peak")})
u3.peak <- apply(refgene.marks.bi["u3", , ], 2, function (x) {Geno(refgene.tx, x, "peak")})

colnames(tx.peak) <- paste(colnames(tx.peak), "tx.peak", sep = ".")
colnames(tss1k.peak) <- paste(colnames(tss1k.peak), "tss1k.peak", sep = ".")
colnames(pm2k.peak) <- paste(colnames(pm2k.peak), "pm2k.peak", sep = ".")
colnames(pm5k.peak) <- paste(colnames(pm5k.peak), "pm5k.peak", sep = ".")
colnames(ex.peak) <- paste(colnames(ex.peak), "ex.peak", sep = ".")
colnames(it.peak) <- paste(colnames(it.peak), "it.peak", sep = ".")
colnames(en.peak) <- paste(colnames(en.peak), "en.peak", sep = ".")
colnames(u5.peak) <- paste(colnames(u5.peak), "u5.peak", sep = ".")
colnames(u3.peak) <- paste(colnames(u3.peak), "u3.peak", sep = ".")

#--- build complete data frame
mel.df <- data.frame(row.names = rownames(refgene.tx), rna.df[-idx, ],
  tx.length, tss1k.length, pm2k.length, pm5k.length, ex.length, it.length, en.length, u5.length, u3.length,
  tx.peak, tss1k.peak, pm2k.peak, pm5k.peak, ex.peak, it.peak, en.peak, u5.peak, u3.peak,
  tx.intensity, ex.intensity, it.intensity, u5.intensity, u3.intensity,
  tx.signal, tss1k.signal, pm2k.signal, pm5k.signal, ex.signal, it.signal, en.signal, u5.signal, u3.signal)

mel.df[is.na(mel.df)] <- 0
save(marks.bi, mel.df, file = "~/Dropbox/Network/R1/melDf.rdt")  # ESB4
