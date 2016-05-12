library(dplyr)
library(GOstats)
library(KEGG.db)
library(org.Hs.eg.db)
library(Category)
library(pathview)
library(ggplot2)

rm(list = ls())
hpc <- "/data/xwang/ADSP"; github <- "~/Dropbox/GitHub/Adsp"

setwd(github); load("Shiny/table.rdt")
setwd(github); load("Shiny/beta/beta.rdt"); geneId <- beta$summary$query

# GO AND KEGG ENRICHMENT 
universe_go <- get("org.Hs.egGO") %>% Lkeys
universe_kegg <- get("org.Hs.egPATH") %>% Lkeys
entrezId <- mget(geneId, org.Hs.egSYMBOL2EG, ifnotfound = NA) %>% na.omit %>% unlist

myGO <- lapply(c("BP", "MF", "CC"), function(category) {
  cat(category, "\n")
  params <- new("GOHyperGParams", geneIds = entrezId, universeGeneIds = universe_go, 
    annotation = "org.Hs.eg.db", ontology = category, pvalueCutoff = 0.001, testDirection = "over")  
  hyperGTest(params) %>% summary
}); names(myGO) <- c("BP", "MF", "CC")

params <- new("KEGGHyperGParams", geneIds=entrezId, universeGeneIds = universe_kegg, 
  annotation="org.Hs.eg.db", categoryName = "KEGG", pvalueCutoff = 0.05, testDirection="over")
kegg <- hyperGTest(params) %>% summary

glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(x) {y <- mget(x, envir=org.Hs.egSYMBOL); paste(y, collapse=";")})
kegg$Symbols <- glist[as.character(kegg$KEGGID)]

# IREGULON 
genes <- unique(c(geneId, alzGene))
genes <- genes[! genes == "-"]
write(genes, file = "Pathway/genes.txt")

# ************************
alzgene <- read.delim("./Alzgene/alzgene20120917.tsv", stringsAsFactors = F)
table(unique(table.gwc$SYMBOL) %in% alzgene$Gene)
table(table.gwc$ID %in% alzgene$Polymorphism)
write(alzgene$Gene, file = "./iRegulon/geneID_Alz.txt")

load("./GWAS/alzMyfit.rdt")
dt <- alzMyfit[c("ID", "PAR", "LOD")]
dt$SIGN <- sign(dt$PAR)
dt$cons <- rep("Protective", nrow(dt))
dt$cons[dt$SIGN == 1] <- "Risky"

pdf("./GWAS/effect_ALZ.pdf")
ggplot(dt, aes(x = LOD, y = abs(PAR))) + 
  geom_point(aes(shape = factor(cons), colour = factor(cons))) +
  theme_bw() + xlab("LOD") + ylab("Effect") +
  theme(panel.border = element_rect(size = 1, color = "grey30")) +
  scale_colour_manual(values = c("chartreuse3", "firebrick1")) +
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold", vjust = 1)) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()
