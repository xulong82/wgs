library(dplyr)
library(ggplot2)
library(ggrepel)

rm(list = ls()) 
setwd("~/Dropbox/GitHub/wgs")

load("./Manu/table1.rdt")
gwas = table1

gwas$Sign <- as.factor(sign(gwas$mean))
gwas$Type <- "NS"
gwas$Type[1:4] <- "S"

gwas$label = ""
gwas$label[1:4] = c("rs74944275", "rs10490263", "rs149372995", "rs140233081")
gwas$label[gwas$label == ""] = NA

gwas$MAF[gwas$MAF > 0.5] = 1 - gwas$MAF[gwas$MAF > 0.5]

pdf("./Manu/top2.pdf", width = 7, height = 4)

ggplot(gwas, aes(x = MAF, y = abs(mean))) + 
  geom_point(data = gwas[1:4, ], aes(size = -log10(P)), shape = 1) + 
# geom_point(aes(colour = Sign, size = -log10(P), shape = Type), alpha = 0.8) +
  geom_point(aes(colour = Sign, size = -log10(P)), alpha = 0.7) +
  scale_size(range = c(3, 8), breaks = seq(6, 9, 0.5)) +
  theme_bw() + xlab("MAF") + ylab("Absolute Effect Size") +
  guides(color = F, text = F, shape = F) + 
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) + 
  theme(axis.text = element_text(size = 15),
        axis.title= element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15),
#       legend.position = "top", 
        legend.key = element_blank()) +
  geom_label_repel(aes(label = label), nudge_x = c(0.1, -0.01, 0.1, 0.1), nudge_y = c(0, 0.2, 0, 0))

dev.off()

gwas_coding$Label <- gwas_coding$Symbol
gwas_coding$Label <- gsub("AC010336", "", gwas_coding$Label)
gwas_coding$Label <- gsub("ZNF684", "", gwas_coding$Label)
gwas_coding$Label[duplicated(gwas_coding$Label)] <- ""
gwas_coding$Sign <- as.factor(sign(gwas_coding$pSnp))

gwas.lod$Sign <- as.factor(sign(gwas.lod$pSnp))
