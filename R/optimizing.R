library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())

arg <- commandArgs(TRUE) 
chr <- paste("chr", arg, sep = "")

load("~/Dropbox/GitHub/glmm/data/mdata.rdt")
load("~/Dropbox/GitHub/glmm/data/kinship.rdt")

load(paste0("/data/xwang/Adsp/golden/geno_", chr, ".rdt"))

load("/data/xwang/glmm/prune.rdt")
prune_v <- prune$V1

load("/data/xwang/glmm/noK_p0.01.rdt")
noK_v <- rownames(noK_p0.01)

geno <- geno[rownames(geno) %in% c(prune_v, noK_v), ]

cov <- mdata[c("Age", "Sex")]
cov$Sex <- as.integer(cov$Sex) - 1
Ad <- as.numeric(mdata$AD2)
g0 <- rep(0, length(Ad))

Sigma <- kinship[[paste0("no_", chr)]]
Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

model <- stan_model("~/Dropbox/GitHub/glmm/Stan/glmm2.stan")
dat = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, Ad = Ad)
dat$sigma = 0.5

# fit <- optimizing(model, data = dat)
# init = list(p = 0, beta = fit$par[paste0("beta[", 1:2, "]")], c = fit$par[paste0("c[", 1:3, "]")])
# init$Sigma = fit$par["Sigma"]
# init$z = fit$par[paste0("z[", 1:576, "]")]
# save(init, file = "~/Dropbox/GitHub/glmm/data/init_glmm.rdt")
load("~/Dropbox/GitHub/glmm/data/init_glmm.rdt")

myGWAS_optimizing <- function(geno) {

  N1 = nrow(geno)
  vId = rownames(geno)
  pId = c("p", paste0("beta[", 1:2, "]"), paste0("c[", 1:3, "]"))
  y1 = matrix(nrow = N1, ncol = 8, dimnames = list(vId, c(pId, "se", "lp")))

  for (i in 1:N1) {
    dat1 = within(dat, {g = geno[i, ]})
    fit1 = optimizing(model, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", init = init, data = dat1)
    se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
    lp1 = sum(fit1$par[paste0("lp[", 1:576, "]")])
    y1[i, ] = c(fit1$par[pId], se1, lp1)
  }
  return(y1)
}

n.core <- 20 
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1 
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
genoList <- mclapply(1:n.core, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

y <- mclapply(genoList, myGWAS_optimizing, mc.cores = n.core)
fit <- do.call(rbind, y)

save(fit, file = paste0("/data/xwang/glmm/glmm2/", chr, ".rdt"))

