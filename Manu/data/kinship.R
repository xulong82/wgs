library(ggplot2)
library(pheatmap)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")
load("./data/kinship.rdt")

pdf("./Pdf/kin-heatmap.pdf", width = 6, height = 4)

pheatmap(cor(kinship$kin23), treeheight_col = 0, treeheight_row = 1e2)

dev.off()

pedigree = c(rep("In", 1545), rep("Out", 164055))
gdt <- data.frame(kinship = kinship$kin23[, "autosome"], pedigree = pedigree)

pdf("./Manu/kinship.pdf", width = 4.5, height = 4, family = "Helvetica")

ggplot(gdt, aes(kinship)) + 
  geom_density(aes(color = pedigree), size = 0.5, show_guide = F) +
  stat_density(aes(kinship, color = pedigree), geom = "line", size = 1, position = "identity") +
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) +
  theme_bw() + xlab("Kinship Coefficient") + ylab("Density") +
  theme(panel.border = element_rect(size = 1, color = "grey30"),
        axis.text = element_text(size = 13),
      	axis.title = element_text(size = 15),
        legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), 
        legend.key = element_blank())

dev.off()

load("./data/KS_doqtl.rdt")
load("./data/KS_emma.rdt")
compare <- cbind(EMMA = c(kin1), DOQTL = c(kin2), KING = c(kinship$autosome))
pheatmap(cor(compare))
