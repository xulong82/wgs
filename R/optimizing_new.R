library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())
setwd("/data/xwang/byglmm")

arg <- commandArgs(TRUE) 
chr <- paste("chr", arg, sep = "")

load("./mdata.rdt")
load("./kinship.rdt")
load("./init_c.rdt")

# load(paste0("../Adsp/golden/geno_", chr, ".rdt"))

load(paste0("./geno/", chr, ".rdt"))
load(paste0("./prior/", chr, ".rdt"))

# load("./prune.rdt")
# geno <- geno[rownames(geno) %in% prune$V1, ]

model <- stan_model("./model_new.stan") # ordered

cov <- mdata[c("Age", "Sex")]
cov$Sex <- as.integer(cov$Sex) - 1

Ad <- as.numeric(mdata$AD2)

g0 <- rep(0, length(Ad))
p0 <- 0

Sigma <- kinship[[paste0("no_", chr)]]
Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

dat = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, prior = p0, Ad = Ad)

myGWAS_optimizing <- function(geno) {

  N1 = nrow(geno)
  vId = rownames(geno)
  p1 = as.matrix(igap[vId, ])
  p1 = p1[, 1] / p1[, 2]

  y1 = matrix(nrow = N1, ncol = 2, dimnames = list(vId, c("p", "se")))

  for (i in 1:N1) {
    dat1 = within(dat, {g = geno[i, ]; prior = p1[i]})
    fit1 = optimizing(model, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", init = init, data = dat1)
    
    se1 = tryCatch(sqrt(diag(solve(-fit1$hessian)))["p"], error=function(e) NULL)
    y1[i, ] = c(fit1$par["p"], se1)
  }
  return(y1)
}

n.core <- 15 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
genoList <- mclapply(1:n.core, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

y <- mclapply(genoList, myGWAS_optimizing, mc.cores = n.core)
fit <- do.call(rbind, y)

save(fit, file = paste0("./c_prior_hessian/", chr, ".rdt"))

