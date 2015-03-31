#--- permutations
Permutation <- function (x) {
  maxLogLik <- 0
  for (i in 1: ncol(pairs)) {
    x1 <- geno[, pairs[, i]]
    #   x1 <- apply(geno[, pairs[, i]], 2, sample)
    data1 <- data.frame(rna = x, x1)
    model1 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "+"))), data1)  # additive
    model2 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "*"))), data1)  # complete
    loglik1 <- logLik(model2) - logLik(model1)
    if (loglik1 >  maxLogLik) maxLogLik <- loglik1
  }
  return(maxLogLik)
}

cat(date(), "--- multiple cores run \n")
n.core = 20
n.loop = 10 
rnaList <- list()
for (i in 1:n.core) rnaList[[i]] <- sample(rna2)
maxLogLik <- rep(0, n.core * n.loop)
for (i in 1:n.loop) {  # n.core by n.loop permutations
  y <- mclapply(rnaList, Permutation, mc.cores = n.core)
  start <- (i - 1) * n.core + 1
  end <- i * n.core 
  maxLogLik[start:end] <- unlist(y)
}
cat(date(), "--- fi \n")

arg <- commandArgs(TRUE)
fname <- paste("permute", arg, sep = "")
assign(fname, maxLogLik)
save(list = fname, file = paste(paste("~/Dropbox/Network/R1/permute", fname, sep = "/"), "rdt", sep = "."))
