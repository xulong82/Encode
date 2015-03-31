# code deleted

#--- archieve --- CODE BELOW IS STEPWISE FORWARD REGRESSION
# model1 <- paste("y ~ ", paste(marks, collapse= "+"))  # additive
# model2 <- paste("y ~ ", paste(marks, collapse= "*"))  # complete
# model3 <- paste(paste("y ~ (", paste(marks, collapse = "+")), ")^2")  # 2nd order
# 
# step.aic <- matrix(nrow = n.epi, ncol = 4, dimnames = list(id.epi, c("mark", "pval", "estimate", "adj.r2")))
# for (i in 1:n.epi) {  # ET1
#   if (i %% 1e3 == 0) cat(paste(as.integer(i * 100 / n.epi), "%", sep = ""), "\n")
#   dat <- data.frame(et = rna.svd$v[1, ], t(epi.geno[, , i]))
# # lmlr <- lm(as.formula(model2), dat) 
#   step1 <- lm(et ~ 1, dat)
#   step2 <- stepAIC(step1, as.formula(model2), direction = "forward", trace = F)
#   coefs <- summary(step2)$coefficients
#   if (max(summary(step2)$df) < 10) {
#     pval <- min(coefs[, "Pr(>|t|)"][-1])
#     mark <- rownames(coefs)[coefs[, "Pr(>|t|)"] == pval]
#     estimate <- coefs[coefs[, "Pr(>|t|)"] == pval, "Estimate"]
#     adj.r2 <- summary(step2)$adj.r.squared
#     step.aic[i, ] = c(mark, pval, estimate, adj.r2)
#   }
# }
# 
# step.aic1 <- step.aic[!is.na(step.aic[, "pval"]), ]
# hist(as.numeric(step.aic1[, "pval"]))
# step.aic2 <- step.aic1[as.numeric(step.aic1[, "pval"]) < 0.01, ]
# table(step.aic2[, "mark"])
#
# #----------------------------------------------------------------------------------------
# myLength <- function(x) {  # Calculate length of consecutive 1 sequences
#   # x: binary sequence; y: length of consecutive 1 sequences
#   y0 <- ifelse(x[1] == 0, 0, 1)  # status: 1
#   id <- ifelse(x[1] == 0, 0, 1)  # id for sequences
#   y1 <- rep(NA, 100000)  # initiate vecotor for length of 1 sequences
#   pb <- txtProgressBar(min = 0, max = 100, style = 3)  # Progress bar
#   len <- length(x)
#   for (i in 2:len) {
#     if (i %% 100000 == 0) {
#       progress <- (i * 100) %/% length(x)
#       setTxtProgressBar(pb, progress)  
#     }  # Update progress bar
#     turn <- x[i] - x[i-1]
#     if (y0) {
#       if (turn == 0) {
#         y0 = y0 + 1
#         if (i == len) y1[id] = y0
#       } else {
#         y1[id] = y0
#         y0 = 0
#       } 
#     } else {
#       if (turn == 1) {
#         y0 = 1
#         id = id + 1
#         if (i == len) y1[id] = y0
#       }
#       y0 <- ifelse(turn == 0, 0, 1)
#     }
#   }
#   close(pb)  # Close progress bar
#   return(y1)
# }
# #----------------------------------------------------------------------------------------
# myCount1 <- function(x) {  # Calculate number of consecutive 1 sequences
#   # x: binary sequence; y: number of consecutive 1 sequences
#   y <- ifelse(x[1] == 0, 0, 1)  # initiate y
#   pb <- txtProgressBar(min = 0, max = 100, style = 3)  # Progress bar
#   for (i in 2:length(x)) {
#     if (i %% 100000 == 0) {
#       progress <- (i * 100) %/% length(x)
#       setTxtProgressBar(pb, progress)  
#     }  # Update progress bar
#     turn <- x[i] - x[i-1]
#     y <- y + ifelse(turn == 1, 1, 0)
#   }
#   close(pb)  # Close progress bar
#   return(y)
# }
# #----------------------------------------------------------------------------------------
# myCount2 <- function(x1, x2) {  # Calculate number of 1 or sum of signal in a region
#   # x1: BED columns 1-3 (1:chrom; 2:chromStart; 3:chromEnd)
#   # x2: Binary genome/epigenome features
#   # y:  number of 1 or sum of signal in x2 with regions defined by x1
#   y <- c(rep(0, length(x1)))
#   x1$chrom <- gsub("chrX", "chr20", x1$chrom)
#   x1$chrom <- gsub("chrY", "chr21", x1$chrom)  # Encode has no Y chromosome
#   x1$chrom <- gsub("X", "20", x1$chrom)
#   x1$chrom <- gsub("Y", "21", x1$chrom)  # Encode has no Y chromosome
#   x1$chrom <- as.numeric(gsub("chr", "", x1$chrom))
#   pb <- txtProgressBar(min = 0, max = 100, style = 3)  # Progress bar  
#   for (i in 1:nrow(x1)) {
#     if (i %% 1000 == 0) {
#       progress <- (i * 100) %/% nrow(x1)
#       setTxtProgressBar(pb, progress) 
#     }  # Update progress bar
#     x1.1 <- x1[i, ]  # Single record
#     offset <- ifelse(x1.1$chrom == 1, 0, sum(chrom.length[1:(x1.1$chrom - 1)]))
#     start <- offset + x1.1$chromStart %/% interval
#     end <- offset + x1.1$chromEnd %/% interval
#     y[i] <- sum(x2[start:end])
#   }
#   close(pb)  # Close progress bar
#   return(y)
# }
# #----------------------------------------------------------------------------------------
# myPattern <- function(x1, x2) {  # Assign pattern of marks to genes
#   # x1: epigenome features
#   # x2: pattern matrix
#   y <- matrix(0, nrow = nrow(x1), ncol = nrow(x2))  # pattern of marks to genes
#   for (i in 1:nrow(x1)) {
#     if (i %% 1000 == 0) print(i)
#     x1.1 <- x1[i, ]  # Single record
#     x1.1[x1.1 > 1] <- 1
#     index <- which(apply(x2, 1, function (x) {paste(x, collapse = "") == paste(x1.1, collapse = "")}))
#     y[i, index] <- 1  # flexible
#   }
#   return(y)
# }
# #----------------------------------------------------------------------------------------

#--- STOP --- CODE BELOW TREATING EACH GENETIC LOCATION A SAMPLE
reg0 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "+"))), data1)  # additive
reg0 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse= "*"))), data1)  # complete

n.start <- 2
n.end <- 6
n.combo <- (n.end - 1) * (n.end - 2) /2
reg1 <- matrix(nrow = n.combo, ncol = 4, dimnames = list(1:n.combo, c("m1", "m2", "lrt", "pval")))
k <- 0
for (i in n.start:(n.end - 1)) {
  for (j in (i + 1): n.end) {
    k = k + 1
    cat(k, paste(colnames(data1)[i], colnames(data1)[j], sep = ":"), "\n")
    dt.tmp <- data1[, c(1, i, j)]
    g0 <- lm(as.formula(paste("rna ~ ", paste(colnames(dt.tmp)[-1], collapse = "+"))), data = dt.tmp)
    g1 <- lm(as.formula(paste(paste("rna ~ (", paste(colnames(dt.tmp)[-1], collapse = "+")), ")^2")), data = dt.tmp)
    #   g0 <- lm(as.formula(paste("rna ~ ", paste(colnames(data1)[-1], collapse = "+"))), data = data1)
    #   g1 <- lm(as.formula(paste(paste("rna ~ ", paste(colnames(data1)[-1], collapse = "+")),
    #                             paste("+", paste(colnames(dt.tmp)[-1], collapse = "*")))), data = data1)
    lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
    pval <- 1 - pchisq(lrt, df = 1)  # Return p value
    reg1[k, ] <- c(i, j, lrt, pval)
  }
}  #  LRT extreme large and p extreme low: R2 very low: eventhough model improves significantly (high LRT, 

#--- RNA ~ marks: stepAIC --------------------------------------------
reg2.step1 <- lm(rna ~ 1, data1) 
reg2.step2 <- stepAIC(reg2.step1, as.formula(paste("~. +", paste(colnames(data1)[-1], collapse = "+"))),
                      direction="forward")  # additive
reg2.step3 <- stepAIC(reg2.step1, as.formula(paste("~. +", paste(colnames(data1)[-1], collapse = "*"))),
                      direction="forward")  # additive & first order interactive
#--- marks ~ marks ---------------------------------------------------
data2 <- log2(marks.1.tx[rowSums(marks.1.tx) > 128, ] + 1)
n.start <- 2
n.end <- ncol(data2) 
n.combo <- (n.end - 1) * (n.end - 2) /2
reg2 <- array(data = NA, c(n.end, n.combo, 4), 
              dimnames = list(colnames(data2), 1:n.combo, c("m1", "m2", "lrt", "pval")))
for (l in 1:ncol(data2)) {
  k <- 0
  dt.tmp1 <- cbind(y = data2[, l], data2[, -l])
  for (i in n.start:(n.end - 1)) {
    for (j in (i + 1): n.end) {
      k = k + 1
      cat(k, paste(colnames(dt.tmp1)[i], colnames(dt.tmp1)[j], sep = ":"), "\n")
      dt.tmp2 <- dt.tmp1[, c(1, i, j)]
      g0 <- lm(as.formula(paste("y ~ ", paste(colnames(dt.tmp2)[-1], collapse = "+"))), data = dt.tmp2)
      g1 <- lm(as.formula(paste(paste("y ~ (", paste(colnames(dt.tmp2)[-1], collapse = "+")), ")^2")), data = dt.tmp2)
      # g0 <- lm(as.formula(paste("y ~ ", paste(colnames(data1)[-1], collapse = "+"))), data = data1)
      # g1 <- lm(as.formula(paste(paste("y ~ ", paste(colnames(data1)[-1], collapse = "+")),
      #                           paste("+", paste(colnames(dt.tmp)[-1], collapse = "*")))), data = data1)
      lrt <- 2 * (logLik(g1) - logLik(g0))  # Likelihood ratio test
      pval <- 1 - pchisq(lrt, df = 1)  # Return p value
      reg2[l, k, ] <- c(i, j, lrt, pval)
    }
  }
}  #  LRT extreme large and p extreme low
#*** A PROBLEM OF LINEAR REGRESSION IN TERMS OF DATA SIZE AND VARIABLE SIZE ***
data3 <- ch12.marks.bi[rowSums(ch12.marks.bi) >= 2, ]

#-- FT: delete Absent patterns
data2 <- cbind(rna = log2(mel.rna$TPM + 1), pt.marks.1.tx)
#-- FT: delete non-expressed, non-marked genes
data2 <- data2[!((data2$V1 == 1) & (data2$rna == 0)), ]  

reg2.step1 <- lm(rna ~ 1, data2)
reg2.step2 <- stepAIC(reg2.step1, as.formula(paste("~. +", paste(colnames(data2)[-1], collapse = "+"))),
                      direction="forward")  # additive
reg2.step3 <- stepAIC(reg2.step1, as.formula(paste(paste("~. + (", paste(colnames(data2)[-1], collapse = "+")), ")^2")),
                      direction="forward")  # additive & first order interaction 
# reg2.step4 <- stepAIC(reg3.step1, as.formula(paste("~. +", paste(colnames(data2)[-1], collapse = "*"))),
#                       direction="forward")  # Error: C stack usage is too close to the limit
idx <- gsub("V", "", names(reg2.step3$coefficients)[-1])

tile.dt2 <- data.frame(value = c(as.matrix(pattern[idx, ])), 
                       idx = rep(paste("P", rownames(pattern[idx, ]), sep = ""), 7), 
                       mark = c(apply(as.matrix(colnames(pattern[idx, ])), 1, function (x) {rep(x, length(idx))})))
pdf(file = "~/Dropbox/Encode/Figures/tile2.pdf", height = 3)
ggplot(tile.dt2, aes(x = idx, y = mark, fill = value)) + geom_tile() +
  scale_fill_gradient(low="white", high="black") +
  theme_bw() + xlab("") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        legend.position = "none")
dev.off()

data3 <- pt.marks.1.tx[, -1]  # Remove non-markded pattern 
reg3.step2 <- list()
for (i in 1:ncol(data3)) {
  print(i)
  reg3.step1 <- lm(as.formula(paste(colnames(data3)[i], "~ 1")), data3)
  reg3.step2[[i]] <- stepAIC(reg3.step1, as.formula(paste(paste("~. + (", paste(colnames(data3)[-i], collapse = "+")), ")^2")),
                             direction="forward")
}

#--- define epigenomic hotspots
idx.mk1 <- matrix(nrow = n.genome, ncol = n.mks)  # new binary data per mark
# epi.mk1 <- list()  # define the epi-geno regions for single mark
for (i in 1:n.mks) {
  cat(i, "/", n.mks, marks[i], "\n")  # progress
  idx <- as.logical(colSums(marks.bi[, i, ]))
  idx.mk1[, i] <- as.numeric(idx)  # new merged binary
  # epi.mk1[[i]] <- myRegion(which(idx))  # epi-geno regions
  # rownames(epi.mk1[[i]]) <- paste(marks[i], rownames(epi.mk1[[i]]), sep = "-")
  # epi.mk1[[i]] <- epi.mk1[[i]][(epi.mk1[[i]]$end - epi.mk1[[i]]$start) >= 1, ]  # length-1
}
dimnames(idx.mk1) <- list(1:n.gnm, marks)
names(epi.mk1) <- marks

#-------
myDel1 <- function (x) {  # delete single 1 in binary sequence
  # x: binary sequence; y: binary sequence w/o single 1
  leng <- length(x)
  y <- x  # return sequence
  if (x[1] - x[2] == 1) y[1] = 0  # head
  if (x[leng] - x[leng-1] == 1) y[leng] = 0  # end
  for (i in 2:(leng - 1)) {
    if (i %% .5e6 == 0) cat(paste(as.integer(i * 100 / leng), "%", sep = ""), "\n")
    if ((x[i] - x[i-1] == 1 ) & (x[i] - x[i+1] == 1)) {
      y[i] = 0
    }
  }
  return(y)
}