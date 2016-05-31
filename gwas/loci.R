library(dplyr)
library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/wgs/loci")

load("2:40973289.rdt")
load("5:102726073.rdt")
load("7:580540.rdt")

sapply(fit, length)

mcmc <- lapply(fit, function(x) { sapply(x, function(y) y["p", ]) })

mcmc <- t(do.call(cbind, mcmc)) %>% as.data.frame
mcmc$P <- pnorm(abs(unlist(mcmc$mean)), sd = unlist(mcmc$sd), lower.tail = F) * 2

mcmc$CHR = gsub(":.*", "", rownames(mcmc))
mcmc$POS = gsub(".*:(.*)_.*", "\\1", rownames(mcmc))

plot(mcmc$POS, -log10(mcmc$P))
