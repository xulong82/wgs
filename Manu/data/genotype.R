library(dplyr)
library(ggplot2)

rm(list = ls()) 
setwd("~/Dropbox/GitHub/Adsp")

# genotype graph

load("./data/mdata.rdt")
load("./data/genoList.rdt")

geno = genoList$glm
geno = genoList$glm_rare

all(colnames(geno) == mdata$SRR)

vars = rownames(geno)
vars = vars[ vars != "" ]

colSums(geno)
plot(colSums(geno), xlab = "Sample Index", ylab = "Alternative allele number")
which(colSums(geno) > 30)
mdata[which(colSums(geno) > 30), ]

sapply(vars, function(x) table(geno[x, ], mdata$AD2))

# risky
var <- "10-120789413" # NANOS1
var <- "19-17830077"
var <- "13-101326992"
var <- "19-7964859"
var <- "10-75148067"
var <- "19-17830077"
var <- "2-191951361"

# protective
var <- "19-4538599" # LGR1
var <- "4-119610606"

dat1 <- mutate(mdata, geno = geno[var, ])
dat1$epi <- with(dat1, paste(Apoe2, geno, sep = ":")) 
dat1$epi <- with(dat1, paste(Apoe4, geno, sep = ":")) 
dat1$epi <- gsub("^", "A", dat1$epi)
dat1$epi <- gsub(":", ":V", dat1$epi)

mycol <- c("grey70", "dodgerblue3", "firebrick1")

pdf(paste0("./geno_graph/", var, ".pdf"), width = 6, height = 2.5, family = "Helvetica")

ggplot(dat1, aes(x = AD2, fill = as.factor(geno))) + 
  geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("")  + coord_flip() +
  scale_fill_manual(values = mycol, breaks=c("0", "1", "2"), labels=c("0/0", "0/1", "1/1")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.key = element_blank()) 

dev.off()

mycol <- c("grey70", "darkorchid2", "dodgerblue3", "chartreuse3", "gold1", "firebrick1")
# mycol <- c("grey70", "grey30", "darkorchid2", "dodgerblue3", "chartreuse3", "gold1", "firebrick1")

pdf(paste0("./geno_graph/", var, "_Apoe4.pdf"), width = 6, height = 2.5, family = "Helvetica")

ggplot(dat1, aes(x = AD2, fill = epi)) + 
  geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("")  + coord_flip() +
  scale_fill_manual(values = mycol) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.key = element_blank()) 

dev.off()
