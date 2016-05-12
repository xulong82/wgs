library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())

chr <- commandArgs(TRUE) 

# setwd("~/Dropbox/GitHub/Adsp/new")
setwd("/data/xwang/adsp3")

load("./mdata.rdt")
load(paste0("/data/xwang/adsp3/R/chr", chr, ".rdt"))
geno <- geno[, mdata$ADSP.Sample.ID]

glm <- stan_model("./glm2.stan")

cov <- mdata[c("Age", "Sex")]
dat = list(N = 570, D = 2, cov = cov)

dat$Ad = rep(0, 570)
dat$Ad[mdata$AD1 %in% c("Probable", "Definite")] = 1
dat$g = rep(0, 570)

# fit <- optimizing(glm, data = dat)
# init_glm <- list(a = fit$par["a"], p = 0, beta = fit$par[paste0("beta[", 1:2, "]")])
# save(init_glm, file = "~/Dropbox/GitHub/Adsp/new/init_glm2.rdt")
load("./init_glm2.rdt")

myGWAS_optimizing <- function(geno) {

  N1 = nrow(geno)
  vId = rownames(geno)
  pId = c("p", "a", paste0("beta[", 1:2, "]"))
  y1 = matrix(nrow = N1, ncol = 5, dimnames = list(vId, c(pId, "se")))

  for (i in 1:N1) {
    dat1 = within(dat, {g = geno[i, ]})
    fit1 = optimizing(glm, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", init = init_glm, data = dat1)
    se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
    y1[i, ] = c(fit1$par[pId], se1)
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

save(fit, file = paste0("/data/xwang/adsp3/glm/chr", chr, ".rdt"))

