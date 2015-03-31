# #--- define epigeno values on hotspots for single mark
# myHPC1 <- function (x) {  # x: hotspots per mark
#   n.epi <- nrow(x)
#   id.epi <- rownames(x)
#   y <- matrix(nrow = n.epi, ncol = n.sps, dimnames = list(id.epi, samples))
#   mark <- gsub("-.*", "", id.epi[1])
#   marks.bi1 <- t(marks.bi[, mark, ])
#
#   for (i in 1:n.sps) {
#     cat(i, "/", n.sps, samples[i], "\n")
#     for (j in 1:n.epi) {
#       if (j %% 1e4 == 0) cat(j, "/", n.epi, "\n")
#       start <- x[j, 1]
#       end <- x[j, 2]
#       y[j, i] <- sum(marks.bi1[start:end, i])
#     }
#   }
#   return(y)
# }
# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network 
# Rev: July 18, 2014

### PEARSON'S CORRELATION TEST

library(parallel)

load("/data/xwang/Network/R/rna.rdt")
load("/data/xwang/Network/R/epigeno.rdt")

arg <- commandArgs(TRUE)
fname <- paste("pearson", arg, sep = "")
  
start <- (as.numeric(arg) - 1) * 4e2 + 1  
if (as.numeric(arg) != 10) {
  end <- (as.numeric(arg)) * 4e2
} else {
  end <- nrow(rna.f)
}

rna <- split(rna.f, rownames(rna.f))  # data frame to list
rna <- rna[start:end]

n.sps <- dim(epi.geno)[1]
n.mks <- dim(epi.geno)[2]
n.epi <- dim(epi.geno)[3]
samples <- dimnames(epi.geno)[[1]]
marks <- dimnames(epi.geno)[[2]]
id.epi <- dimnames(epi.geno)[[3]]

myCor <- function (x) {  # Pearson's correlation test
  rna1 <- c(as.matrix(x))  
  y <- array(dim = c(n.mks, n.epi, 2), 
             dimnames = list(marks, id.epi, c("pho", "p")))
  pho1 <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, marks))
  p1 <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, marks))
  for (j in 1:n.mks) {  # per mark
    epi1 <- epi.geno[, j, ]
    for (k in 1:n.epi) {  # per epi-spot
      if (k %% 1e4 == 0) cat(k, "in", n.epi, "\n")
      t.cor = cor.test(rna1, epi1[, k], method = "pearson")
      pho1[k, j] = t.cor$estimate
      p1[k, j] = t.cor$p.value
    }
  }
  y[, , 1] <- t(pho1) 
  y[, , 2] <- t(p1) 
  return(y)
}

#-------
time.start <- strptime(date(), "%a %b %d %H:%M:%S %Y")
assign(fname, mclapply(rna, myCor, mc.cores = 20))
save(list = fname, 
     file = paste(paste("/data/xwang/Network/R/scan1", fname, sep = "/"), "rdt", sep = "."))
time.end <- strptime(date(), "%a %b %d %H:%M:%S %Y")
difftime(time.end, time.start, units = "auto")

# #------- Linear regression
# myHPC2 <- function (x) {  # epigeno2
#   epi2 <- epi.geno2  # no filter
#   n.epi <- dim(epi2)[3]
#   id.epi <- dimnames(epi2)[[3]]
#   lod <- matrix(nrow = n.epi, ncol =nrow(x), dimnames = list(id.epi, rownames(x)))
#   for (i in 1:nrow(x)) {
#     cat(i, "in", nrow(x), "genes", "\n")
#     rna.dat <- c(as.matrix(x[i, ]))
#     loglik0 <- logLik(lm(rna.dat ~ 1))
#     loglik1 <- rep(NA, n.epi)
#     for (j in 1:n.epi) {
#       if (j %% 1e4 == 0) cat(j, "in", n.epi, "epi", "\n")
#       dat <- data.frame(y = rna.dat, mark = epi2[, , j])
#       loglik1[j] <- logLik(lm(y ~ ., data = dat))
#     }
#     lod[, i] <- (loglik1 - loglik0) / log(10)  # lod
#   }
#   return(lod)
# }
# 
# #-------
# arg <- commandArgs(TRUE)
# fname <- paste("lod2", arg, sep = "-")
# 
# offset <- (as.numeric(arg) - 1) * 1000  # 1000 genes per node
# 
# rna.list <- list()
# for (i in 1:20) {  # core number: 50 genes per core 
#   start = offset + (i - 1) * 50 + 1
#   end = offset + i * 50
#   rna.list[[i]] <- rna.f[start:end, ]
# }
# 
# #-------
# time.start <- strptime(date(), "%a %b %d %H:%M:%S %Y")
# assign(fname, mclapply(rna.list, myHPC2, mc.cores = 20))
# save(list = fname, 
#      file = paste(paste("/data/xwang/Network/R", fname, sep = "/"), "rdt", sep = "."))
# time.end <- strptime(date(), "%a %b %d %H:%M:%S %Y")
# difftime(time.end, time.start, units = "auto")

#-------
# myHPC1 <- function (x) {  # epigeno1
#   lod <- list()
#   for (i in 1:n.mks) {  # per mark
#     cat(i, "in", n.mks, "marks", "\n")
#     epi1 = epi.geno1[[i]]  
#     n.epi = nrow(epi1)
#     id.epi = rownames(epi1)
#     lod1 <- matrix(nrow = n.epi, ncol =nrow(x), dimnames = list(id.epi, rownames(x)))
#     for (j in 1:nrow(x)) {  # per gene
#       cat(j, "in", nrow(x), "genes", "\n")
#       rna.dat <- c(as.matrix(x[j, ]))
#       for (k in 1:n.epi) {  # per epi
#         if (k %% 1e5 == 0) cat(k, "in", n.epi, "epi", "\n")
#         dat <- data.frame(y = rna.dat, mark = epi1[k, ])
#         g0 <- lm(y ~ 1, data = dat)
#         g1 <- lm(y ~ ., data = dat)
#         lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
#         lod1[k, j] <- lrt / (2 * log(10))  # Lod
#       }
#     }
#     lod[[i]] <- lod1
#   }
#   return(lod)
# }
# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network 
# Rev: July 18, 2014

### MARKER AND GENE SELECTION BASED ON EFFECT SIZE (PHO) AND P VALUE

cutoff.pho <- .9
cutoff.pval  <- .01
idx.epi.per.gene <- list()  # post-filter epi-hotspots idx per gene

#--- per gene: epi-spots with required p and pho on 1 or more marks
arg <- commandArgs(TRUE)
name1 <- paste("pearson", arg, sep = "")
name2 <- paste("index", arg, sep = "")
load(paste(paste("/data/xwang/Network/R/scan1", name1, sep = "/"), "rdt", sep = "."))    
data1 <- get(name1)
rm(list = name1)  # release the memory

for (i in 1:length(data1)) {  # genes per rdt file
  if (i %% 1e2 == 0) cat(i, "in", length(data1), "\n")
  temp1 <- data1[[i]]
  temp.p <- temp1[, , "p"]  # p value
  temp.pho <- temp1[, , "pho"]  # pho value
  null <- which(apply(temp.p, 2, function (x) {length(which(is.na(x))) == 4}))
  temp.p <- temp.p[, -null]  # zero variation in epigenotype
  temp.pho <- temp.pho[, -null]
  log.p <- apply(temp.p, 2, function (x) {min(x, na.rm = T) < cutoff.pval})
  log.pho <- apply(temp.pho, 2, function (x) {max(abs(x), na.rm = T) > cutoff.pho})
  log <- log.p & log.pho
  idx.epi.per.gene[[i]] <- unname(which(log))
}
assign(name2, idx.epi.per.gene)
save(list = name2,
     file = paste(paste("/data/xwang/Network/R/scan1", name2, sep = "/"), "rdt", sep = "."))

###########################################################
# idx.epi.per.gene <- list()
# data <- paste("index", 1:10, sep = "")
# pop1 <- 0  # index
# for (k in 1:length(data)) {  # per rdt file
#   cat(k, "in 10 rdt files", "\n")
#   load(paste(paste("/data/xwang/Network/R/scan1", data[k], sep = "/"), "rdt", sep = "."))    
#   data1 <- get(data[k])
#   rm(list = data[k])  # release the memory
#   pop1 <- c(pop1, length(data1))
#   pop2 <- cumsum(pop1[1:k])[k]
#   for (i in 1:length(data1)) {
#     idx.epi.per.gene[[pop2+i]] <- data1[[i]]
#   }
# }
# save(idx.epi.per.gene, file = "/data/xwang/Network/R/scan1/idx.epi.per.gene.rdt")
  
# cutoff.beta <- .5
# epi.per.gene2 <- list()  # beta value per gene
# idx.epi.per.gene2 <- list()  # epi-hotspots ids per gene
# #--- per gene: epi-hotspots with significant beta
# name.beta <- paste("beta", 1:10, sep = "-")
# for (k in 1:length(name.beta)) {  # per rdt file: 1e3 genes
#   cat(k, "in 10 rdt files", "\n")
#   idx.offset <- (k - 1) * 1e3 
# 
#   load(paste(paste("/data/xwang/Network/R/first.scan", name.beta[k], sep = "/"), "rdt", sep = "."))    
#   beta1 <- get(name.beta[k])
#   rm(list = name.beta[k])  # release the memory
# 
#   #--- per gene, identify hotspots with abs(beta) > .5 on 1 or more marks
#   for (i in 1:1e3) {  # 1e3 genes per rdt file
#     temp <- beta1[[i]]
#     idx <- idx.offset + i
#     logit <- apply(t(temp), 2, function (x) {max(abs(x), na.rm = T) > cutoff.beta})
#     epi.per.gene2[[idx]] <- temp[logit, ]
#     idx.epi.per.gene2[[idx]] <- rownames(epi.per.gene2[[idx]])
#   }
# }
# save(epi.per.gene2, file = "/data/xwang/Network/R/epi.per.gene2.rdt")
# save(idx.epi.per.gene2, file = "/data/xwang/Network/R/idx.epi.per.gene2.rdt")
# 
# library(parallel)
# 
# #--- regression: choose epi-spots with good effect sizes
# load("/data/xwang/Network/R/rna.rdt")
# load("/data/xwang/Network/R/epigeno2.rdt")
# load("/data/xwang/Network/R/idx.epi.per.gene.rdt")
# 
# # idx.epi.per.gene <- idx.epi.per.gene[1:1e3]
# 
# # choose 1e4 genes with the top sd
# rna.f <- rna4[apply(rna4, 1, function (x) {max(x) > 3e1}), ]
# rna.f <- log2(rna.f + 1)  # log2 transformation
# rna.f <- rna.f[apply(rna.f, 1, function (x) {max(x) - min(x) > 2}), ]
# rna.f <- rna.f[names(sort(apply(rna.f, 1, sd), decreasing = T))[1:1e4], ]  
# 
# list1 <- list()
# for (i in 1:length(idx.epi.per.gene)) {
#   if (i %% 1e3 == 0) cat(i, "in", length(idx.epi.per.gene), "\n")
#   list2 <- list()
#   list2$rna <- rna.f[i, ]
#   list2$epi <- idx.epi.per.gene[[i]]
#   list1[[i]] <- list2
# }
# 
# myLm1 <- function (x) {  # per gene
#   dt.rna <- c(as.matrix(x$rna))
#   epi <- epi.geno2[, , x$epi]
#   n.mks <- dim(epi)[[2]]
#   n.epi <- dim(epi)[[3]]
#   marks <- dimnames(epi)[[2]]
#   id.epi <- dimnames(epi)[[3]]
#   beta <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, marks))
#   for (i in 1:n.mks) {
#     cat(i, "in", n.mks, "marks", "\n")
#     epi1 <- t(epi[, i, ])
#     for (j in 1:n.epi) {
#       dt.mark <- c(as.matrix(epi1[j, ]))
#       lm1 = lm(dt.rna ~ dt.mark)
#       beta[j, i] <- lm1$coefficients["dt.mark"]
#     }
#   }
#   return(beta)
# }
# 
# arg <- commandArgs(TRUE)
# fname <- paste("beta", arg, sep = "-")
# offset <- (as.numeric(arg) - 1) * 1e3
# list2 <- list1[(offset + 1):(offset + 1e3)]
# 
# time.start <- strptime(date(), "%a %b %d %H:%M:%S %Y")
# assign(fname, mclapply(list2, myLm1, mc.cores = 20))
# time.end <- strptime(date(), "%a %b %d %H:%M:%S %Y")
# difftime(time.end, time.start, units = "auto")
# save(list = fname, file = paste(paste("/data/xwang/Network/R", fname, sep = "/"), "rdt", sep = "."))
# 
