library(dplyr)
library(ggplot2)

rm(list = ls()) 
setwd("~/Dropbox/GitHub/Adsp")

load("./Manu/epis.rdt")
for(obj in names(epis)) assign(obj, epis[[obj]])

vep = epis2$vep
vep = epis4$vep

gwas.vep <- cbind(vep, gwas[match(vep$UID, gwas$UID), ])
sort(table(gwas.vep$Consequence))
gwas.vep[gwas.vep$Consequence == "missense_variant", ]

gwas = sex$gwas
gwas = epis2$gwas
gwas = epis4$gwas

gwas$pSex = gwas$"beta[2]"
gwas$pSnp = gwas$"beta[5]"
gwas$pEpis = gwas$"beta[6]"

x = gwas[c("pApoe2", "pSnp", "pEpis", "LOD2")]  
x$Sign <- as.factor(sign(x$pSnp + x$pApoe2 + x$pEpis))

x = gwas[c("pApoe4", "pSnp", "pEpis", "LOD2")]  
x$Sign <- as.factor(sign(x$pSnp + x$pApoe4 + x$pEpis))

x = gwas[c("pSex", "pSnp", "pEpis", "LOD2")]  
x$Sign <- as.factor(sign(x$pSnp + x$pSex + x$pEpis))

dt = data.frame(melt(x[1:2]), x[3:5])

pdf("./Manu/sex.pdf", width = 6, height = 3, family = "Helvetica")
pdf("./Manu/epis2.pdf", width = 6, height = 3, family = "Helvetica")
pdf("./Manu/epis4.pdf", width = 6, height = 3, family = "Helvetica")

ggplot(dt, aes(x = value, y = pEpis)) + 
  geom_point(aes(colour = Sign, size = LOD2), alpha = 0.9, shape = 111) + 
  facet_grid(. ~ variable) + scale_size(range = c(2, 10)) +
  theme_bw() + xlab("Main Effect") + ylab("Epistasis Effect") + 
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) + 
  guides(color = F, size = F, text = F)

dev.off()

gwas_coding$Label <- gwas_coding$Symbol
gwas_coding$Label <- gsub("AC010336", "", gwas_coding$Label)
gwas_coding$Label <- gsub("ZNF684", "", gwas_coding$Label)
gwas_coding$Label[duplicated(gwas_coding$Label)] <- ""
gwas_coding$Sign <- as.factor(sign(gwas_coding$pSnp))
