# Copyright: Xulong Wang (xulong.wang@jax.org)
# Rev: May 30, 2014

#--- genomic features
chrom.info <- read.delim("~/Dropbox/X/chromInfo.txt", header = F)  # chromosome length
chrom.unit <- 200  # unit: nucleotide length
chrom.leng <- chrom.info$V2 %/% chrom.unit + 1  # chromosome length in chrom.unit
chrom.bina <- rep(0, sum(chrom.leng))

#----------------------------------------------------------------------------------------
Colname <- function (x) {  # Assign colnames for ENCODE broadPeak file
  colnames(x)[1:3] <- c("chrom", "chromStart", "chromEnd") 
  return(x)
}
#----------------------------------------------------------------------------------------
Binary <- function(x) {  # Generate binary file (y) for ENCODE broadPeak file (x: BED column 1-3)
  y <- chrom.bina  # Initialize the return genome
  x$chrom <- gsub("chrX", "chr20", x$chrom)
  x$chrom <- gsub("chrY", "chr21", x$chrom)
  x$chrom <- as.numeric(gsub("chr", "", x$chrom))
  for (i in 1:nrow(x)) {
    if (i %% 1e3 == 0) cat(i, "in", nrow(x), "\n")
    x1 <- x[i, ]  # Single record
    offset <- ifelse(x1$chrom == 1, 0, sum(chrom.leng[1:(x1$chrom - 1)]))
    start <- offset + x1$chromStart %/% chrom.unit
    end <- offset + x1$chromEnd %/% chrom.unit
    y[start:end] <- 1
  }
  return(y)
}
#----------------------------------------------------------------------------------------
Signal <- function(x) {  # Generate peak intensity file (y) from ENCODE broadPeak file
  y <- chrom.bina  # Initialize the return genome
  x$chrom <- gsub("chrX", "chr20", x$chrom)
  x$chrom <- gsub("chrY", "chr21", x$chrom)  # Encode has no Y chromosome
  x$chrom <- as.numeric(gsub("chr", "", x$chrom))
  for (i in 1:nrow(x)) {
    if (i %% 1e3 == 0) cat(i, "in", nrow(x), "\n")
    x1 <- x[i, ]  # Single record
    offset <- ifelse(x1$chrom == 1, 0, sum(chrom.leng[1:(x1$chrom - 1)]))
    start <- offset + x1$chromStart %/% chrom.unit
    end <- offset + x1$chromEnd %/% chrom.unit
    y[start:end] <- x1$V4
  }
  return(y)
}
#----------------------------------------------------------------------------------------
Peak <- function(x) {  # compute number of consecutive 1 (peak) in a binary data
  y <- ifelse(x[1] == 0, 0, 1)  # initiate y
  if (length(x) > 1) {
    for (i in 2:length(x)) {
      turn <- x[i] - x[i-1]
      y <- y + ifelse(turn == 1, 1, 0)
    }
  }
  return(y)
}
#----------------------------------------------------------------------------------------
Geno <- function(x1, x2, x3) {  # compute epigenotypic features 
  # compute x3 (length|peak|signal) in a binary genomic data (x2) in regions defined by bed data (x1)
  y <- c(rep(0, nrow(x1)))
  x1$chrom <- as.character(x1$chrom)
  x1$chrom <- gsub("chrX", "chr20", x1$chrom)
  x1$chrom <- gsub("chrY", "chr21", x1$chrom)
  x1$chrom <- as.numeric(gsub("chr", "", x1$chrom))
  myFun <- function (x) {
    if (x3 == "length") return(sum(x))  # number of 1: total coverage
    if (x3 == "peak") return(Peak(x))  # number of separate peaks
    if (x3 == "signal") return(mean(x[x!=0]))  # mean of the peak intensity
  }
  for (i in 1:nrow(x1)) {
    if (i %% 1e4 == 0) cat(i, "in", nrow(x1), "\n")
    x1.1 <- x1[i, ]  # single record
    offset <- cumsum(c(0, chrom.leng))[x1.1$chrom]
    start <- offset + x1.1$chromStart %/% chrom.unit
    end <- offset + x1.1$chromEnd %/% chrom.unit
    y[i] <- myFun(x2[start:end])
  }
  return(y)
}

#--- Genetic database
refgene.ex <- Colname(read.delim("~/Dropbox/X/refGeneExon.bed", header = F, stringsAsFactors = F))
refgene.it <- Colname(read.delim("~/Dropbox/X/refGeneIntron.bed", header = F, stringsAsFactors = F))
refgene.u5 <- Colname(read.delim("~/Dropbox/X/refGene5pUTR.bed", header = F, stringsAsFactors = F))
refgene.u3 <- Colname(read.delim("~/Dropbox/X/refGene3pUTR.bed", header = F, stringsAsFactors = F))
refgene.tx <- Colname(read.delim("~/Dropbox/X/refGeneTX.bed", header = T, stringsAsFactors = F)[, c(3, 5, 6, 12)])

tx.chrom = aggregate(chrom ~ name2, refgene.tx, function(x) {x[1]})
tx.start = aggregate(chromStart ~ name2, refgene.tx, min)
tx.end = aggregate(chromEnd ~ name2, refgene.tx, max)
refgene.tx <- data.frame(row.names = tx.chrom$name2,  # extended version 
                         chrom = tx.chrom$chrom, chromStart = tx.start$chromStart, chromEnd = tx.end$chromEnd)

refgene.en <- refgene.pm5k <- refgene.pm2k <- refgene.tss1k <- refgene.tx
refgene.tss1k$chromStart <- refgene.tx$chromStart - 1e3  # TSS: -1k -> +1k TSS
refgene.tss1k$chromEnd <- refgene.tx$chromStart + 1e3
refgene.pm2k$chromStart <- refgene.tx$chromStart - 2e3  # Promoter: -2k -> TSS
refgene.pm2k$chromEnd <- refgene.tx$chromStart
refgene.pm5k$chromStart <- refgene.tx$chromStart - 5e3  # Promoter: -5k -> TSS
refgene.pm5k$chromEnd <- refgene.tx$chromStart
refgene.en$chromStart <- refgene.tx$chromStart - 2e5  # Enhancer: -2e2k -> 2e2k TSS
refgene.en$chromEnd <- refgene.tx$chromStart + 2e5

refgene.tx.bi <- Binary(refgene.tx)
refgene.tss1k.bi <- Binary(refgene.tss1k)
refgene.pm2k.bi <- Binary(refgene.pm2k)
refgene.pm5k.bi <- Binary(refgene.pm5k)
refgene.ex.bi <- Binary(refgene.ex)
refgene.it.bi <- Binary(refgene.it)
refgene.en.bi <- Binary(refgene.en)
refgene.u5.bi <- Binary(refgene.u5)
refgene.u3.bi <- Binary(refgene.u3)

cap.ex.bi <- c(0, abs(diff(refgene.ex.bi))) | c(abs(diff(refgene.ex.bi)), 0) | refgene.ex.bi
refgene.en.bi <- as.numeric(refgene.en.bi - refgene.tss1k.bi == 1)  # 1k away from TSS
refgene.en.bi <- as.numeric(refgene.en.bi - as.numeric(cap.ex.bi) == 1)  # 200bp away from exon

refgene.bi <- cbind(refgene.tx.bi, refgene.tss1k.bi, refgene.pm2k.bi, refgene.pm5k.bi, 
  refgene.ex.bi, refgene.it.bi, refgene.en.bi, refgene.u5.bi, refgene.u3.bi)
rm(refgene.tx.bi, refgene.tss1k.bi, refgene.pm2k.bi, refgene.pm5k.bi, 
  refgene.ex.bi, refgene.it.bi, refgene.en.bi, refgene.u5.bi, refgene.u3.bi, cap.ex.bi)
colnames(refgene.bi) <- gsub("refgene.", "", colnames(refgene.bi))
colnames(refgene.bi) <- gsub("\\.bi", "", colnames(refgene.bi))
rownames(refgene.bi) <- 1:nrow(refgene.bi)

save.image("~/Dropbox/Network/R1/encodex.rdt")
