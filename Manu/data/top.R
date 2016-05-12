library(dplyr)
library(ggplot2)

rm(list = ls()) 
setwd("~/Dropbox/GitHub/Adsp")

#

load("./data/glmList.rdt")
for(obj in names(glmList)) assign(obj, glmList[[obj]])

gwas.lod <- filter(gwas, LOD > 15) # permutation cut
gwas.lod$Sign <- as.factor(sign(gwas.lod$pSnp))

# 

load("./Manu/epis.rdt")
for(obj in names(epis)) assign(obj, epis[[obj]])

vep = epis2$vep
gwas.lod = epis2$gwas

vep = epis4$vep
gwas.lod = epis4$gwas

gwas.vep <- cbind(vep, gwas.lod[match(vep$UID, gwas.lod$UID), ])
sort(table(gwas.vep$Consequence))
gwas.vep[gwas.vep$Consequence == "missense_variant", ]

gwas.lod$Sign <- as.factor(sign(gwas.lod$pEpis))

pdf("./Manu/glm.pdf", width = 6, height = 4, family = "Helvetica")

ggplot(gwas.lod, aes(x = MAF, y = pSnp)) + 
  geom_point(aes(colour = Sign, size = LOD), shape = 111) + 
  scale_size(range = c(2, 15)) +
  theme_bw() + xlab("MAF") + ylab("Effect") + xlim(c(0, 0.53)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) + 
  guides(color = F, text = F) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key = element_blank()) 

dev.off()

gwas_coding$Label <- gwas_coding$Symbol
gwas_coding$Label <- gsub("AC010336", "", gwas_coding$Label)
gwas_coding$Label <- gsub("ZNF684", "", gwas_coding$Label)
gwas_coding$Label[duplicated(gwas_coding$Label)] <- ""
gwas_coding$Sign <- as.factor(sign(gwas_coding$pSnp))

gwas.lod$Sign <- as.factor(sign(gwas.lod$pSnp))
