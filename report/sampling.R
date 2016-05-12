library(dplyr)
library(ggplot2)
library(Biobase)

rm(list = ls())
hpc <- "/data/xwang/Adsp"
github <- "~/Dropbox/GitHub/Adsp"

setwd(github)

load("./Manu/sampling.rdt")

vars = sampling$var
name = paste0("c", gsub("-", "p", vars))

setwd(hpc)

load("./GWAS/optimizing.rdt")
gwas <- do.call(rbind, glm.fit)

load("./sampling/c1p567242.rdt")

waic = sapply(c1p567242, function(x) x$waic) %>% unlist
summary = do.call(c, sapply(c1p567242, function(x) x$summary))

opt = gwas[match(names(waic), gwas$UID), ]
rownames(opt) = NULL

sampling$opt = opt
sampling$waic = waic
sampling$summary = summary

setwd(github)
save(sampling, file = "./Manu/sampling.rdt")

