library(dplyr)
library(KEGGREST)
library(biomaRt)
library(ggplot2)

rm(list = ls())
options(stringsAsFactors = F)

alzgenes = list() # output
save(alzgenes, file = "./Manu/alzgenes.rdt")

load("~/Desktop/GWAS/optimizing.rdt")
fit = do.call(rbind, glm.fit)
load("~/Desktop/GWAS/QTLRel.rdt")
lmm = do.call(rbind, lmm.fit)

setwd("~/Dropbox/GitHub/Adsp")
alz <- read.delim("./docs/alzgene20120917.tsv")
var <- gsub("-Women", "", alz$Polymorphism)[2:20]
var <- intersect(var, fit$ID)

fit <- fit[fit$ID %in% var, ]
rownames(fit) = NULL

setwd("/data/xwang/Adsp")

geno <- lapply(1:22, function(chr) {
  cat(chr, "\n")
  load(paste0("./golden/geno_chr", chr, ".rdt"))
  geno[rownames(geno) %in% var, ]
})

geno <- do.call(rbind, geno)

# Post

load("./Manu/alzgenes.rdt")
for(obj in names(alzgenes)) assign(obj, alzgenes[[obj]])

summary(fit$LOD)
colnames(fit)[21] = "pSNP"

all(fit$UID == lmm$UID)
plot(fit$LOD, lmm$LOD, xlim = c(0, 2.5), ylim = c(0, 2.5))
abline(0, 1)

dt = fit
dt$GLMM = dt$LOD
dt$LMM = lmm$LOD

dt$Sign = "Risky"
dt$Sign[dt$pSNP < 0] = "Protective"

pdf("./ALZ/models.pdf", width = 6, height = 4, family = "Helvetica")

ggplot(dt, aes(x = GLMM, y = LMM)) + 
  geom_point(aes(size = MAF, colour = Sign), shape = 111) + 
  scale_size(range = c(2, 15)) +
  theme_bw() + xlab("LOD-GLMM") + ylab("LOD-LMM") + 
  scale_color_manual(values = c("dodgerblue3", "firebrick1")) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(size = .5),
        axis.text = element_text(size = 12),
        axis.title= element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key = element_blank()) 

dev.off()

qplot(MAF, LOD, data = fit, size = abs(pSNP))

xx1 = fit[fit$LOD > 2, ]
alz[alz$Polymorphism %in% xx1$ID, ]
xx1$MAF # Comparable to GMAF

xx2 = mutate(mdata, geno = geno[xx1$UID[3], ])
qplot(x = AD2, data = xx2, geom = "bar", fill = as.factor(geno), position = "fill")
table(xx2$geno, xx2$AD2)

xx3 = as.data.frame.matrix(table(xx2$geno, xx2$AD2))

(WT = xx3[1, ] / colSums(xx3))
(HET = xx3[2, ] / colSums(xx3))
(HOM = xx3[3, ] / colSums(xx3))

(dt = rbind(WT = WT, HET = HET, HOM = HOM))
(dt = melt(as.matrix(dt)))
dt$X1 = factor(dt$X1, levels = c("WT", "HET", "HOM"))
dt$X2 = factor(dt$X2, levels = c("No", "Possible", "Probable", "Definite"))

mycol <- c("grey70", "dodgerblue3", "firebrick1")

myplot <- function(var) {
  cat(var, "\n")
  xx = mutate(mdata, geno = geno[var, ])
  ggplot(xx, aes(x = AD2, fill = as.factor(geno))) + 
    geom_bar(position = "fill", width = 0.65) + ggtitle(var) +
    theme_bw() + xlab("") + ylab("Frequency")  + coord_flip() +
    scale_fill_manual(values = mycol, 
                      breaks=c("0", "1", "2"), 
                      labels=c("0/0", "0/1", "1/1")) +
    theme(panel.border = element_blank(),
          axis.line = element_line(color = 'grey30'),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_blank(), 
          legend.key = element_blank()) 
}

lapply(fit$UID, function(var) {
  filename = paste0("./ALZ/geno/", var, ".pdf")
  pdf(filename, width = 6, height = 2.5, family = "Helvetica")
  print(myplot(var))
  dev.off()
})

pdf("./ALZ/geno.pdf", width = 6, height = 2.5, family = "Helvetica")
lapply(fit$UID, function(var) myplot(var))
dev.off()

qplot(X2, value, data = dt, geom = "line", group = X1, color = X1)

lm(AD1 ~ geno, data = xx2) %>% summary

xx4 = rbind(xx2, xx2)
lm(AD1 ~ geno, data = xx4) %>% summary

# Clear trend, mild LOD, small sample size was a reason

for (i in 1:nrow(xx1)) {
  print(xx1[i, ])
  xx2 = mutate(mdata, geno = geno[xx1$UID[i], ])
  print(table(xx2$geno, xx2$AD2))
}

# Graph

load("./data/hg19.rdt")
source("./GWAS/loci200k.R")

load("~/Desktop/GWAS/optimizing.rdt")
gwas = do.call(rbind, glm.fit)

for (var in fit$UID) {

cat(var, "\n")

filename = paste0("./ALZ/", var, ".pdf")
pdf(filename, width = 9, height = 6)
make.locus.200k("marker", var, gwas, hg19)
dev.off()

}

# GWAS CATELOG

library(gwascat)
(cur = makeCurrentGwascat()) 

gwas <- read.delim("./docs/gwas_catalog_v1.0-downloaded_2015-07-22.tsv")

table(gwas$DISEASE.TRAIT)

gwas_genes <- gwas$MAPPED_GENE %>% unique
gwas_genes <- lapply(gwas_genes, function(x) unlist(strsplit(x, ", "))) 
gwas_genes <- lapply(gwas_genes, function(x) unlist(strsplit(x, " - "))) 
gwas_genes <- unlist(gwas_genes) %>% unique

data$gwas <- gwas

(x = intersect(genes, gwas_genes))
y = lapply(gwas$MAPPED_GENE, function(x) unlist(strsplit(x, ", "))) 
y = lapply(y, function(i) unlist(strsplit(i, " - "))) 
idx = sapply(y, function(i) any(i %in% x))

gwas_select = gwas[idx, ]
table(gwas_select$CONTEXT)
sort(table(gwas_select$DISEASE.TRAIT), decreasing = T)[1:20]

disease = "Rheumatoid arthritis"
disease = "Alzheimer's disease"
disease = "IgG glycosylation"
disease = "Obesity-related traits"
disease = "Type 2 diabetes"

gwas_select %>% filter(grepl(disease, DISEASE.TRAIT))
gwas_select$MAPPED_GENE[grepl(disease, gwas_select$DISEASE.TRAIT)]

gene = "NANOS1"

ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
getBM("name_1006", "external_gene_name", gene, ensembl)

keggFind("pathway", "Long-term potentiation")  # KEGG: LTR

nanos1 = read.csv("NANOS1/metatargetome.csv", stringsAsFactors = F)
nanos1_gk = myGK(nanos1$Target.Gene)
