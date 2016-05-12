library(dplyr)
library(Biobase)

chr = paste0("chr", 1:22)
setwd("/data/xwang/Adsp/permutation")

P300A = lapply(chr, function(x) {
  load(paste0("./R1/", x, ".rdt")); get(x)
}) 

P300A = do.call(cbind, P300A) %>% rowMax

P300B = do.call(cbind, lapply(chr, function(x) {
  load(paste0("./R2/", x, ".rdt")); get(x)
})) %>% rowMax

P300C = do.call(cbind, lapply(chr, function(x) {
  load(paste0("./R3/", x, ".rdt")); get(x)
})) %>% rowMax

permutation = list(P300A, P300B, P300C)

# Graph

library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")

load("./Manu/data.rdt")
load("./data/permutation.rdt")

lp0 = sum(data$par0[paste("lp[", 1:576, "]", sep = "")])
permutation = data$permutation

x = do.call(cbind, permutation)
x = 2 * (x - lp0)

quantiles = c(0.7, 0.8, 0.9, 0.95, 0.99)

a =sapply(quantiles, function(y) quantile(x[, 1], y))
b =sapply(quantiles, function(y) quantile(rowMax(x[, 1:2]), y))
c =sapply(quantiles, function(y) quantile(rowMax(x), y))

dt = rbind(a, b, c)
rownames(dt) = c("2 million", "4 million", "6 million")
dt = melt(dt)
dt$Number = dt$Var1

pdf("./Manu/permutation.pdf", width = 6, height = 4)

ggplot(dt, aes(x = Var2, y = value, group = Number)) + 
  geom_line(aes(color = Number), size = 1) + geom_point(size = 3) +
  theme_bw() + xlab("") + ylab("LOD") +
  scale_color_manual(values = c("dodgerblue3", "firebrick1", "chartreuse3")) + 
  theme(panel.border = element_rect(size = 1, color = "grey30"),
        legend.key = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16))

dev.off()

y = rowMax(x)

summary(y)
quantile(y, 0.95)

init = data$model$init 
init = colMeans(init)

permutation = data$model$permutation

x1 = sapply(permutation, function(x) quantile(x, 0.95))
x2 = x1 - init

plot(x)

y1 = sapply(1:22, function(x) 2 * (permutation[[x]] - init[x]))
quantile(c(y1), 0.90)

dt = data.frame(LOD = c(y1))

pdf("./Pdf/permutation.pdf", width = 5, height = 4, family = "Helvetica")

ggplot(dt, aes(x = LOD)) +
  geom_density(color = "dodgerblue3", size = 1) +
  geom_vline(x = 19, linetype = 2, size = 1, color = "firebrick1") +
  theme_bw() + xlab("") + ylab("Density") +
  theme(panel.border = element_rect(size = 1, color = "grey30"),
        axis.text = element_text(size = 13),
      	axis.title = element_text(size = 15, vjust = 1))

dev.off()
  
