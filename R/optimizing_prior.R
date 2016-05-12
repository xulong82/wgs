library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())

arg <- commandArgs(TRUE) 
chr <- paste("chr", arg, sep = "")

load("~/Dropbox/GitHub/glmm/data/mdata.rdt")
load("~/Dropbox/GitHub/glmm/data/kinship.rdt")

load(paste0("/data/xwang/glmm/geno/", chr, ".rdt"))
load(paste0("/data/xwang/glmm/prior/", chr, ".rdt"))

load("/data/xwang/glmm/prune.rdt")
prune_v <- prune$V1

igap_p <- pnorm(abs(igap$Beta), sd = igap$SE, lower.tail = F) * 2
igap_v <- rownames(igap)[igap_p < 0.01]

load("/data/xwang/glmm/noK_p0.01.rdt")
noK_v <- rownames(noK_p0.01)

geno <- geno[rownames(geno) %in% c(prune_v, igap_v, noK_v), ]

cov <- mdata[c("Age", "Sex")]
cov$Sex <- as.integer(cov$Sex) - 1
Ad <- as.numeric(mdata$AD2)
g0 <- rep(0, length(Ad))

Sigma <- kinship[[paste0("no_", chr)]]
Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

# model <- stan_model("~/Dropbox/GitHub/glmm/Stan/glmm_prior.stan")
model <- stan_model("~/Dropbox/GitHub/glmm/Stan/glmm_prior2.stan")
dat = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, prior = 0, Ad = Ad)

# fit <- optimizing(model, data = dat)
# init = list(beta = fit$par[paste0("beta[", 1:2, "]")], c = fit$par[paste0("c[", 1:3, "]")])
# init$sigma = fit$par["sigma"]; init$t = init$p = 0
# init$sigma = fit$par["sigma"]; init$sigma2 = fit$par["sigma2"]; init$t = init$p = 0
# init$z = fit$par[paste0("z[", 1:576, "]")]
# save(init, file = "~/Dropbox/GitHub/glmm/data/init_glmm_prior.rdt")
load("~/Dropbox/GitHub/glmm/data/init_glmm_prior.rdt")

myGWAS_optimizing <- function(geno) {

  N1 = nrow(geno)
  vId = rownames(geno)
  pId = c("p", paste0("beta[", 1:2, "]"), paste0("c[", 1:3, "]"))
  y1 = matrix(nrow = N1, ncol = 7, dimnames = list(vId, c(pId, "se")))
  p1 = igap[vId, 1] / igap[vId, 2]

  for (i in 1:N1) {
    dat1 = within(dat, { g = geno[i, ]; prior = p1[i] })
    fit1 = optimizing(model, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", init = init, data = dat1)
    se1 = tryCatch(sqrt(diag(solve(-fit1$hessian)))["p"], error=function(e) NULL)
    y1[i, ] = c(fit1$par[pId], se1)
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

save(fit, file = paste0("/data/xwang/glmm/glmm_prior2/", chr, ".rdt"))

