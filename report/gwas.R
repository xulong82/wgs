library(dplyr)
library(ggplot2)
library(Biobase)

rm(list = ls())
hpc <- "/data/xwang/Adsp"
github <- "~/Dropbox/GitHub/Adsp"

setwd(github)

load("./data/genInitList.rdt")
loglik0 <- sum(par0[paste("logLik[", 1:576, "]", sep = "")])

load("./data/sampling.rdt")
loglik0 <- sum(par0[paste("logLik[", 1:576, "]", sep = "")])

load("./Manu/model.rdt")
init <- model$init %>% colMeans

chr <- paste0("chr", 1:22)
chr <- paste0("chr", c(1:22, "X", "Y"))

names(init) = chr

setwd(hpc)

# GLM

glm.fit <- lapply(chr, function(chr) {
  cat(chr, "\n")
# load(paste("./GLM/biVar/", chr, ".rdt", sep = ""))
# load(paste("./GLM/biVar_rare/", chr, ".rdt", sep = ""))
# load(paste("./GLM/epis_Apoe2/", chr, ".rdt", sep = ""))
# load(paste("./GLM/epis_Apoe4/", chr, ".rdt", sep = ""))
  load(paste("./epistasis/sex/", chr, ".rdt", sep = ""))
# load(paste("./Meta/meta_", chr, ".rdt", sep = ""))
# load(paste0("./optimizing/", chr, ".rdt"))
# load(paste0("./null/", chr, ".rdt"))
  load(paste0("./golden/meta_", chr, ".rdt"))
  fit <- get(chr)
# fit <- cbind(meta[match(names(fit), meta$UID), ], lp = fit)
  fit <- cbind(meta[match(rownames(fit), meta$UID), ], fit)
  fit$LOD <- 2 * (fit$lp - init[chr])
# fit$LOD <- 2 * (fit$loglik - init[chr])
  rownames(fit) <- NULL; fit
}); names(glm.fit) <- chr

# QTLRel

lmm.fit <- lapply(chr, function(chr) { cat(chr, "\n")
# load(paste("./QTLRel/", chr, ".rdt", sep = ""))
  load(paste("./QTLRel_rare/", chr, ".rdt", sep = ""))
  load(paste("./Meta/meta_", chr, ".rdt", sep = ""))
  fit.p <- lapply(get(chr), function (x) x$p) %>% unlist
  fit.v <- lapply(get(chr), function (x) x$v) %>% unlist
  fit.e <- lapply(get(chr), function (x) do.call(rbind, x$parameters))
  fit.e <- do.call(rbind, fit.e)
  fit <- cbind(meta[match(names(fit.p), meta$UID), ], fit.e)
  mutate(fit, LOD = fit.p, V = fit.v)
}); names(lmm.fit) <- chr

# LMM

hlm.fit <- lapply(chr, function(chr) { cat(chr, "\n")
  load(paste("./LMM/", chr, ".rdt", sep = ""))
  load(paste("./Meta/meta_", chr, ".rdt", sep = ""))
  fit <- get(chr)
  cbind(meta[match(rownames(fit), meta$UID), ], fit)
}); names(hlm.fit) <- chr

save(glm.fit, file = "./GWAS/null.rdt")
save(glm.fit, file = "./GWAS/optimizing.rdt")
save(glm.fit, file = "./GWAS/GLM.rdt")
save(glm.fit, file = "./GWAS/GLM_rare.rdt")
save(lmm.fit, file = "./GWAS/QtlRel.rdt")
save(lmm.fit, file = "./GWAS/QtlRel_rare.rdt")
save(hlm.fit, file = "./GWAS/Hlm.rdt")
save(glm.fit, file = "./GWAS/GLM_epis_Apoe2.rdt")
save(glm.fit, file = "./GWAS/GLM_epis_Apoe4.rdt")
save(glm.fit, file = "./GWAS/GLM_epis_Sex.rdt")

# Permutation

glm.permutation <- lapply(chr, function(chr) { cat(chr, "\n")
  load(paste("./GLM/permutation/", chr, ".rdt", sep = ""))
  get(chr)
}); names(glm.permutation) <- chr

cut <- quantile(unlist(glm.permutation), 0.90)

quantile(gwas$LOD, 0.999)
gwas <- gwas[gwas$LOD > quantile(gwas$LOD, 0.999), ] # 0.001 quantile

# MODEL COMPARISON

# MAF > 0.01

load("GWAS/optimizing"); glm <- do.call(rbind, glm.fit)
load("GWAS/QtlRel.rdt"); lmm <- do.call(rbind, lmm.fit)
load("GWAS/Hlm.rdt"); hlm <- do.call(rbind, hlm.fit)
load("GWAS/GLM_epis_Apoe2.rdt"); glm_epis2 <- do.call(rbind, glm.fit)
load("GWAS/GLM_epis_Apoe4.rdt"); glm_epis4 <- do.call(rbind, glm.fit)

lodALL <- cbind(GLMM = glm$LOD, LMM = lmm$LOD, LMER = hlm$LMER_LOD, CLMM = hlm$CLMM_LOD)
lodALL <- as.data.frame(lodALL)
lodALL$epi_Apoe2 <- glm_epis2$LOD
lodALL$epi_Apoe4 <- glm_epis4$LOD
save(lodALL, file = "./GWAS/lodALL.rdt")

# MAF < 0.01

load("GWAS/GLM_rare.rdt"); glm <- do.call(rbind, glm.fit)
load("GWAS/QtlRel_rare.rdt"); lmm <- do.call(rbind, lmm.fit)

random <- sample(1:2e7, 10)
all(glm$MAF[random] == lmm$MAF[random])

lodALL <- cbind(GLM = glm$LOD, LMM = lmm$LOD)
save(lodALL, file = "./GWAS/lodALL(rare).rdt")

# Quantiles

load("./GWAS/lodALL.rdt")
lodALL = lod[c("GLMM", "LMM")] 
lodALL = lod[c("GLMM", "LMM", "LMER", "CLMM")]

quantiles = c(0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995)

lod = lapply(quantiles, function(y) {
  index = apply(lodALL, 2, function(x) x > quantile(x, y))
  lodALL[rowSums(index) > 0, ]
})

names(lod) = paste0("q", quantiles)

sapply(lod, nrow)
lapply(lod, cor)

setwd(github)
save(lod, file = "./Manu/lod.rdt")

