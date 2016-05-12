library(rstan)
library(parallel)
library(dplyr)

rm(list = ls())
setwd("~/Dropbox/GitHub/glmm")

load("./data/mdata.rdt") # phenotype
load("./data/kinship.rdt") # kinship
load("./Apoe/apoe.rdt") # genotype

cov <- mdata[c("Age", "Sex")]
cov$Sex <- as.integer(cov$Sex) - 1
Ad <- as.numeric(mdata$AD2)
Sigma <- kinship[["no_chr19"]]
Sigma[Sigma < 0] <- 0
L <- t(chol(Sigma))

load("./Apoe/prior.rdt")
g0 <- rep(0, length(Ad))
geno <- g1 <- apoe[rownames(apoe) %in% rownames(igap), ]
p1 <- as.matrix(igap[rownames(g1), ])
N1 <- nrow(g1)
vId <- rownames(g1)
pos <- gsub("*-", "", vId) %>% as.numeric

n.core <- 8
u.core <- round(nrow(geno) / n.core)
idx1 <- ((1:n.core) - 1) * u.core + 1
idx2 <- c((1:(n.core - 1)) * u.core, nrow(geno))
genoList <- mclapply(1:n.core, function(x) geno[idx1[x]:idx2[x], ], mc.preschedule = F)

# GLMM with prior

model <- stan_model("./Stan/glmm_prior.stan") 
data = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, prior = 0, Ad = Ad)

zc = optimizing(model, data = data)
init = list(p = 0, t = 0, beta = zc$par[paste0("beta[", 1:2, "]")], c = zc$par[paste0("c[", 1:3, "]")], z = zc$par[paste0("z[", 1:576, "]")])
init$sigma = 0.25

my_optimizing = function(g1) {
  N1 = nrow(g1); vId = rownames(g1)
  p1 <- p1[vId, 1] / p1[vId, 2]
  pId <- c("p", "t", "sigma")
  fit <- matrix(nrow = N1, ncol = 4, dimnames = list(vId, c(pId, "se")))
  for (i in 1:N1) {
    dat1 = within(dat, { g = g1[i, ]; prior = p1[i] })
    fit1 = optimizing(model, hessian = T, algorithm = "LBFGS", init = init, data = dat1)
    se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
    fit[i, ] = c(fit1$par[pId], se1)
  }
return(fit) }

y <- mclapply(genoList, myopt, mc.cores = n.core)
fit <- do.call(rbind, y)

save(fit, file = "./Apoe/myfit.rdt")

fit = as.data.frame(fit)
x4 = -log10(pnorm(abs(fit$p), sd = fit$se, lower.tail = F) * 2)

pval = pnorm(abs(fit$mu), sd = fit$sigma, lower.tail = F) * 2
pval = pnorm(abs(fit$p), sd = fit$se, lower.tail = F) * 2

load("./gwas/c_wo_prior.rdt")
wo_prior <- z[rownames(fit)]
x3 = -log10(pchisq(wo_prior * 2, df = 1, lower.tail = F))

all(rownames(p1) == rownames(fit))
x1 = -log10(pnorm(abs(p1[, 1]), sd = p1[, 2], lower.tail = F) * 2)
data = data.frame(prior_p = x1, prior_t = p1[, 1] / p1[, 2], wprior = x2, wo_prior_lrt = x3, wo_prior = x4)
plot(data$wo_prior_lrt, data$wo_prior)

plot(pos, data$wprior)
plot(pos, data$wo_prior)

data2 = data[order(data$wprior, decreasing = T), ]
head(data2, n = 20)

x = c("19-45411941", "19-45412079")
data[x, ]

load("./Apoe/snp2.rdt")

all(mdata$SRR == colnames(y$geno))
mdata = cbind(mdata, t(y$geno))

dat1 <- mutate(mdata, geno = y$geno[2, ])
ggplot(dat1, aes(x = AD2, fill = as.factor(geno))) + geom_bar(position = "fill", width = 0.65) 

# Stan for linear model method in PLINK

model <- stan_model("./Stan/lm.stan") 
dat = list(N = 576, D = 2, cov = cov,  g = g0, Ad = mdata$AD1)

zc = optimizing(model, data = dat, hessian = T)
zc$par
sqrt(diag(solve(-zc$hessian)))

init = list(a = 0.078, p = 0, beta = zc$par[paste0("beta[", 1:2, "]")], sigma = 0.25)

pId <- c("a", "p", "sigma")
fit <- matrix(nrow = N1, ncol = 4, dimnames = list(vId, c(pId, "se")))
  
for (i in 1:N1) {
  dat1 = within(dat, { g = g1[i, ] })
  fit1 = optimizing(model, hessian = T, algorithm = "LBFGS", init = init, data = dat1)
  se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
  fit[i, ] = c(fit1$par[pId], se1)
}

fit <- as.data.frame(fit)
pval <- pnorm(abs(fit$p), sd = fit$se, lower.tail = F) * 2
plot(pos, -log10(pval))
abline(h = 7.3, col = "red")

# Actual PLINK output
# Stan ouputs on LM are the same as PLINK and SNPTEST

load("./Apoe/plink.rdt")
plot(plink2$BP, -log10(plink2$P))
abline(h = 7.3, col = "red")

# Stan for categorical model 
# LM to GLM will drop significance, why and how much?

model <- stan_model("./Stan/noK.stan") 
data = list(N = 576, K = 4, D = 2, cov = cov,  g = g0, Ad = Ad)

zc = optimizing(model, data = data, hessian = T); zc$par
sqrt(diag(solve(-zc$hessian)))

init = list(p = 0, beta = zc$par[paste0("beta[", 1:2, "]")], c = zc$par[paste0("c[", 1:3, "]")])

pId <- c("p", "beta[1]", "beta[2]", "c[1]", "c[2]", "c[3]")
fit <- matrix(nrow = N1, ncol = 7, dimnames = list(vId, c(pId, "se")))

for (i in 1:N1) {
  dat1 = within(data, { g = g1[i, ] })
  fit1 = optimizing(model, hessian = T, algorithm = "LBFGS", init = init, data = dat1)
  se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
  fit[i, ] = c(fit1$par[pId], se1)
}

glm.fit <- fit <- as.data.frame(fit)
pval <- pnorm(abs(fit$p), sd = fit$se, lower.tail = F) * 2
plot(pos, -log10(pval), ylim = c(0, 7.5))
abline(h = 7.3, col = "red")

# LM to LMM in Stan
# GLM to GLMM drops significance, why?

model <- stan_model("./Stan/glmm.stan") 
data = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, Ad = Ad)

zc = optimizing(model, hessian = T, data = data); zc$par
sqrt(diag(solve(-zc$hessian)))
init = list(p = 0, beta = zc$par[paste0("beta[", 1:2, "]")], c = zc$par[paste0("c[", 1:3, "]")], z = zc$par[paste0("z[", 1:576, "]")])

my_optimizing = function(g1) {
  N1 <- nrow(g1); vId = rownames(g1);
  pId <- c("p", "beta[1]", "beta[2]", "c[1]", "c[2]", "c[3]")
  fit <- matrix(nrow = N1, ncol = 7, dimnames = list(vId, c(pId, "se")))
  for (i in 1:N1) {
    dat1 = within(data, { g = g1[i, ] })
    fit1 = optimizing(model, hessian = T, algorithm = "LBFGS", init = init, data = dat1)
    se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
    fit[i, ] = c(fit1$par[pId], se1)
  }
return(fit) }

y <- mclapply(genoList, my_optimizing, mc.cores = n.core)
fit <- do.call(rbind, y)

glmm.fit <- fit <- as.data.frame(fit)
glmm_test.fit <- fit <- as.data.frame(fit)
pval <- pnorm(abs(fit$p), sd = fit$se, lower.tail = F) * 2
plot(pos, -log10(pval), ylim = c(0, 7.5))
abline(h = 7.3, col = "red")

# P-values drop by p or se(p)? If so, why?

test = list(glm = glm.fit, glmm = glmm.fit)
save(test, file = "./Apoe/test.rdt")

summary(glmm.fit$se - glm.fit$se)
summary(glmm.fit$p - glm.fit$p)

glm.fit$P = -log10(pnorm(abs(glm.fit$p), sd = glm.fit$se, lower.tail = F) * 2)
glmm.fit$P = -log10(pnorm(abs(glmm.fit$p), sd = glmm.fit$se, lower.tail = F) * 2)

summary(glm.fit$P - glmm.fit$P)
index = which.max(glm.fit$P - glmm.fit$P)

glm.fit[index, ]
glmm.fit[index, ]

g1[index, ]
dat1 <- mutate(mdata, geno = g1[index, ])
ggplot(dat1, aes(x = AD2, fill = as.factor(geno))) + geom_bar(position = "fill", width = 0.65) 

# what was the ramdom effect

model <- stan_model("./Stan/glmm_test2.stan") 

data = list(N = 576, K = 4, D = 2, cov = cov, L = L, g = g0, Ad = Ad)
fit0 = optimizing(model, hessian = T, data = data)
se0 = sqrt(diag(solve(-fit0$hessian)))["p"]
-log10(pnorm(fit0$par["p"], sd = se0, lower.tail = F) * 2)

sam0 = sampling(model, data = data, chains = 3, cores = 3, iter = 600, warmup = 200)

u0 = fit0$par[paste0("u[", 1:576, "]")]
cor(mdata$AD1, u0)

# high correlation suggests necessity of adding the random term
# the random term can affect estimating variant effect in two directions, obserbing the effect, and increase resolution

L2 = 0.25 * L # this prior has a meaning
# random term models association between gentoype and phenotype
# predict the phenotypic similarity between individuals
# use the covariance structure to implement

dat1 = within(data, { g = g1[index, ] })
fit1 = optimizing(model, hessian = T, data = dat1)
se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
-log10(pnorm(fit1$par["p"], sd = se1, lower.tail = F) * 2)

Sr = genPositiveDefMat("eigen",dim = 576, lambdaLow = 0.1)$Sigma
Lr <- t(chol(Sr))

data = list(N = 576, K = 4, D = 2, cov = cov, L = Lr, g = g1[index, ], Ad = Ad)
fit1 = optimizing(model, hessian = T, data = data)
se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
-log10(pnorm(fit1$par["p"], sd = se1, lower.tail = F) * 2)

z1 = fit1$par[paste0("u[", 1:576, "]")]
cor(mdata$AD1, z1)
