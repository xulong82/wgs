library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())

chr <- commandArgs(TRUE) 

load("~/Dropbox/GitHub/Adsp/new/mdata.rdt")
load(paste0("/data/xwang/adsp3/R/chr", chr, ".rdt"))
geno <- geno[, mdata$ADSP.Sample.ID]

glm <- stan_model("~/Dropbox/GitHub/Adsp/new/glm.stan")

cov <- mdata[c("Age", "Sex")]
dat = list(N = 570, K = 4, D = 2, cov = cov)
dat = within(dat, { Ad = as.numeric(mdata$AD1); g = rep(0, 570) })

# fit <- optimizing(glm, data = dat)
# init_glm <- list(p = 0, beta = fit$par[paste0("beta[", 1:2, "]")], c = fit$par[paste0("c[", 1:3, "]")])
# save(init_glm, file = "~/Dropbox/GitHub/Adsp/new/init_glm.rdt")
load("~/Dropbox/GitHub/Adsp/new/init_glm.rdt")

myGWAS_optimizing <- function(geno) {

  N1 = nrow(geno)
  vId = rownames(geno)
  pId = c("p", paste0("beta[", 1:2, "]"), paste0("c[", 1:3, "]"))
  y1 = matrix(nrow = N1, ncol = 8, dimnames = list(vId, c(pId, "se", "lp")))

  for (i in 1:N1) {
    dat1 = within(dat, {g = geno[i, ]})
    fit1 = optimizing(glm, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", init = init_glm, data = dat1)
    se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
    lp1 = sum(fit1$par[paste0("lp[", 1:570, "]")])
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

save(fit, file = paste0("/data/xwang/adsp3/glm/chr", chr, ".rdt"))

