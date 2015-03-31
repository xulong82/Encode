# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Epigenomic Network 
# Rev: August 18, 2014

library(parallel)

load("/data/xwang/Network/R/rna.mod.rdt")
load("/data/xwang/Network/R/epigeno.rdt")

cat(date(), "--- PEARSON'S CORRELATION TEST \n")

rna.mod <- rna.mod[, c(2:11, 1)]  # match the samples
rna <- split(rna.mod, 1:nrow(rna.mod))  # df to list

n.mks <- dim(epigeno)[2]
n.epi <- dim(epigeno)[3]
id.mks <- dimnames(epigeno)[[2]]
id.epi <- dimnames(epigeno)[[3]]

myAns <- function (x) {  # Pearson's correlation test
  y <- array(dim = c(n.mks, n.epi, 5), 
             dimnames = list(id.mks, id.epi, c("pho", "p", "q", "beta", "logLik")))
  pho1 <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, id.mks))
  p1 <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, id.mks))
  q1 <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, id.mks))
  beta1 <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, id.mks))
  logLik1 <- matrix(nrow = n.epi, ncol = n.mks, dimnames = list(id.epi, id.mks))
  model1 <- lm(x ~ 1)
  for (j in 1:n.mks) {  # per mark
    cat(j, "in", n.mks, "\n")
    epi1 <- t(epigeno[, j, ])
    for (k in 1:n.epi) {  # per epi-spot
      t.cor = cor.test(x, epi1[k, ], method = "pearson")
      model2 = lm(x ~ epi1[k, ])
      pho1[k, j] = t.cor$estimate
      p1[k, j] = t.cor$p.value
      beta1[k, j] = model2$coefficients[2]
      logLik1[k, j] = logLik(model2) - logLik(model1)
    }
    q1[, j] <- p.adjust(p1[, j], method = "fdr")
  }
  y[, , 1] <- t(pho1) 
  y[, , 2] <- t(p1) 
  y[, , 3] <- t(q1) 
  y[, , 4] <- t(beta1) 
  y[, , 5] <- t(logLik1) 
  return(y)
}

#-------
scan1 <- mclapply(rna, myAns, mc.cores = 20)
save(scan1, file = "/data/xwang/Network/R/scan1/scan1.rdt")

cat(date(), "--- MARKER SELECTION BASED ON FDR(Q) VLUE \n")

cutoff.qval <- .05
scan1f <- list()
for (i in 1:length(scan1)) {
  if (i %% 1e2 == 0) cat(i, "in", length(scan1), "\n")
  scan11 <- scan1[[i]]
  null <- which(apply(scan11[, , "q"], 2, function (x) {length(which(is.na(x))) == 4}))
  scan11 <- scan11[, -null, ]
  idx <- apply(scan11[, , "q"], 2, function (x) {min(x, na.rm = T) < cutoff.qval})
  scan1f[[i]] <- scan11[, idx, ]
}
save(scan1f, file = "/data/xwang/Network/R/scan1/scan1f.rdt")

# cat(date(), "--- MARKER AND GENE SELECTION BASED ON EFFECT SIZE (PHO) AND P VALUE \n")
# cutoff.pho <- .5
# cutoff.pval  <- .05
# idx.epi <- list()
# 
# #--- per gene: epi-spots with required p and pho on 1 or more marks
# load("/data/xwang/Network/R/scan1/scan1.rdt")
# for (i in 1:length(sscan1)) {
#   if (i %% 1e2 == 0) cat(i, "in", length(sscan1), "\n")
#   temp1 <- sscan1[[i]]
#   temp.p <- temp1[, , "p"]  # p value
#   temp.pho <- temp1[, , "pho"]  # pho value
#   temp.beta <- temp1[, , "beta"]  # beta
#   temp.logLik <- temp1[, , "logLik"]  # logLik
# 
#   null <- which(apply(temp.p, 2, function (x) {length(which(is.na(x))) == 4}))
#   if (length(null) != 0) {
#     temp.p <- temp.p[, -null]  # zero-variation in epigenotype
#     temp.pho <- temp.pho[, -null]
#     temp.beta <- temp.beta[, -null]
#     temp.logLik <- temp.logLik[, -null]
#   }
#   #--- significant spots in any mark
#   true.p <- apply(temp.p, 2, function (x) {min(x, na.rm = T) < cutoff.pval})
#   true.pho <- apply(temp.pho, 2, function (x) {max(abs(x), na.rm = T) > cutoff.pho})
#   idx.epi$any[[i]] <- unname(which(true.p & true.pho))
#   scan1.f[[i]] <- scan1[, , which(true.p & true.pho)] 
#   #--- significant spots per mark
#   true.p <- apply(temp.p, 1, function(x) {x < cutoff.pval})
#   true.pho <- apply(temp.pho, 1, function(x) {abs(x) > cutoff.pho})
#   idx.epi$pmark[[i]] <- apply(true.p & true.pho, 2, function (x) {unname(which(x))})
#   #--- spots with 2 or more significant marks
#   idx.mmark <- 0
#   for (j in 1:ncol(temp.p)) {
#     temp.p1 <- sort(temp.p[, j])
#     temp.pho1 <- sort(abs(temp.pho[, j]), decreasing = T)
#     if (length(temp.p1) > 1) {
#       if ((temp.p1[2] < cutoff.pval) & (temp.pho1[2] > cutoff.pho)) 
#         idx.mmark <- c(idx.mmark, j)
#     }
#   }
#   idx.epi$mmark[[i]] <- idx.mmark[-1]
# }
# cat(date(), "fi \n")
# save(idx.epi, file = "/data/xwang/Network/R/scan1/idx.epi.rdt")

