library(dplyr)
library(ggplot2)

rm(list = ls())
setwd("/data/xwang/Adsp")

# QQ: Null vs Additive
load("./GWAS/null.rdt")
null = sort(do.call(rbind, glm.fit)$LOD)
load("./GWAS/optimizing.rdt")
real = sort(do.call(rbind, glm.fit)$LOD)

# Epistasis

load("./GWAS/optimizing.rdt")
glm0 = do.call(rbind, glm.fit)
load("./GWAS/GLM_epis_Sex.rdt")
glmSex = do.call(rbind, glm.fit)

all(glm0$UID == glmSex$UID) 
glmSex$LOD2 = 2 * (glmSex$lp - glm0$lp)

gwas = glm0[glm0$LOD > 15, ]
gwas = glmSex[glmSex$LOD2 > 15, ]

load("./GWAS/GLM.rdt")
glm0 = do.call(rbind, glm.fit)
load("./GWAS/GLM_epis_Apoe2.rdt")
glm2 = do.call(rbind, glm.fit)
load("./GWAS/GLM_epis_Apoe4.rdt")
glm4 = do.call(rbind, glm.fit)

glm0 = glm0[glm0$CHR %in% 1:22, ]

all(glm0$UID == glm2$UID) 
all(glm0$UID == glm4$UID) 

glm2$LOD2 = 2 * (glm2$loglik - glm0$loglik)
glm4$LOD2 = 2 * (glm4$loglik - glm0$loglik)

gwas = glm2[glm2$LOD2 > 15, ]
gwas = glm4[glm4$LOD2 > 15, ]

load("./VEP/vepList.rdt")

vep <- do.call(rbind, lapply(vepList, function (x) x[x$X.Uploaded_variation %in% gwas$UID, ]))
vep$Symbol <- gsub(".*SYMBOL=(.*)", "\\1", vep$Symbol)
rownames(vep) <- rownames(gwas) <- NULL
colnames(vep)[grep("Uploaded", colnames(vep))] <- "UID"

epis = list()

epis$epis2$vep = vep
epis$epis2$gwas = gwas

epis$epis4$vep = vep
epis$epis4$gwas = gwas

setwd("~/Dropbox/GitHub/Adsp")
save(epis, file = "./Manu/epis.rdt")

null = sort(glm0$LOD)
real = sort(glm2$LOD)
real = sort(glm4$LOD)

dt = data.frame(null, real)

setwd("~/Dropbox/GitHub/Adsp")

png("./Manu/qq.png", width = 1200, height = 1000, res = 200)
png("./Manu/qq2.png", width = 1200, height = 1000, res = 200)
png("./Manu/qq4.png", width = 1200, height = 1000, res = 200)

ggplot(dt, aes(x = null, y = real)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw() + xlab("w/o interaction") + ylab("w interaction") +
# theme_bw() + xlab("NULL") + ylab("FULL") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "firebrick1") +
  xlim(c(0, 26)) + ylim(c(0, 26)) +
  theme(panel.border = element_rect(size = 1, color = "grey30"),
  axis.text = element_text(size = 13),
  axis.title = element_text(size = 16))

dev.off()

