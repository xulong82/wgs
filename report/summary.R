library(dplyr)
library(ggplot2)
library(pheatmap)
library(biomaRt)
library(stargazer)

rm(list = ls()) 
setwd("~/Dropbox/GitHub/Adsp")

load("Manu/data.rdt")
for(obj in names(data)) assign(obj, data[[obj]])

load("./Manu/epis.rdt")

load("data/glmList.rdt")
load("data/glmList_rare.rdt")
load("data/glmList_epis_Apoe2.rdt")
load("data/glmList_epis_Apoe4.rdt")

load("data/optimizing.rdt")
load("data/lmmList.rdt")

list <- glmList
list <- lmmList
list <- optimizing

for(obj in names(list)) assign(obj, list[[obj]])

gwas_lod <- filter(gwas, LOD > 15) # permutation cut
gwas_lod <- filter(gwas, LOD > 19) # permutation cut
gwas_lod <- filter(gwas, LOD > quantile(gwas$LOD, 0.99)) # 0.01%

vep_lod <- filter(vep, UID %in% gwas_lod$UID)

lm(gwas_lod$LOD.LMM ~ gwas_lod$LOD)
plot(gwas_lod$LOD, gwas_lod$LOD.LMM, xlab = "GLM (Stan)", ylab = "LMM (QTLRel)")
abline(a = 3.211, b = 1.239)

group <- list()
group[[1]] <- gwas_lod[1, ]
group_idx <- 1

for (i in 2:nrow(gwas_lod)) {
  chromosome = gwas_lod$CHR[i] == gwas_lod$CHR[i-1]
  position = gwas_lod$POS[i] - gwas_lod$POS[i-1] < 1e6
  other = abs(gwas_lod$pSnp[i] - gwas_lod$pSnp[i-1]) < 1
  
  if ( all(chromosome, position, other) )
    group[[group_idx]] = rbind(group[[group_idx]], gwas_lod[i, ])
  else {
    group_idx = group_idx + 1
    group[[group_idx]] = gwas_lod[i, ]
  }
}

group <- lapply(group, function(x) cbind(PEAK = x$UID[which.max(x$LOD)], x))
group <- do.call(rbind, group)

var <- as.character(unique(group$PEAK))
var <- var[! grepl("X", var)]

gwas_vep_lod <- cbind(vep_lod, groupAll[match(vep_lod$UID, groupAll$UID), ])
gwas_vep_lod <- cbind(vep_lod, gwas_lod[match(vep_lod$UID, gwas_lod$UID), ])
gwas_vep_lod <- gwas_vep_lod[with(gwas_vep_lod, pSnp + pApoe4 + pEpis > 1), ]
gwas_vep_lod <- gwas_vep_lod[with(gwas_vep_lod, pSnp + pApoe2 + pEpis < -1), ]

varAll <- gwas$UID[! grepl("[X|Y]", gwas$UID)]
geno_select <- read.plink("/data/xwang/ADSP/Plink/autosome.bed", select.snps = varAll)

geno_select_lod <- which(varAll %in% gwas.lod$UID)
LD_select <- ld(geno_select$genotypes[, geno_select_lod], stats = c("D.prime", "R.squared"), depth = 100)
LD.DP <- as.matrix(LD_select$D.prime)
LD.R2 <- as.matrix(LD_select$R.squared)
diag(LD.R2) <- diag(LD.DP) <- 1
LD.DP[lower.tri(LD.DP)] <- t(LD.DP)[lower.tri(LD.DP)]
LD.R2[lower.tri(LD.R2)] <- t(LD.R2)[lower.tri(LD.R2)]

var_select <- rownames(LD.DP)
gwas.lod_ld <- gwas.lod[match(var_select, gwas.lod$UID), ]
chr <- gwas.lod_ld$CHR; pos <- gwas.lod_ld$POS
idx = sapply( 1:(length(var_select)-1), function (i) {
  sapply( (i + 1) : length(var_select), function (j) {
    if ( LD.DP[i, j] > 0.5 & chr[i] == chr[j] & abs(pos[i] - pos[j]) < 1e6 )
      ifelse( gwas.lod_ld$LOD[i] > gwas.lod_ld$LOD[j], j, i)
})}) %>% unlist %>% unique
gwas.lod_ld <- gwas.lod_ld[-idx, ]

# protein coding variants

so_coding <- names(cons)[1:6]
vep_coding <- vep[apply(sapply(so_coding, function(x) grepl(x, vep$Consequence)), 1, any), ]
vep_coding <- vep_coding[! duplicated(paste0(vep_coding$UID, vep_coding$Symbol)), ]
gwas_coding <- cbind(vep_coding, gwas[match(vep_coding$UID, gwas$UID), ])
gwas_coding <- gwas_coding[order(gwas_coding$LOD, decreasing = T), ]

gwas_coding <- gwas_coding[with(gwas_coding, abs(pSnp) > 0.2), ]
gwas_coding <- gwas_coding[with(gwas_coding, pSnp + pApoe2 + pEpis < -1), ]
gwas_coding <- gwas_coding[with(gwas_coding, pSnp + pApoe4 + pEpis > 1), ]
rownames(gwas_coding) <- NULL

geneId <- vep_coding$Symbol %>% unique
summary <- queryMany(geneId, scopes="symbol", species="human", fields = c("name", "summary"))
summary <- summary[! duplicated(summary$query), c("query", "name", "summary")] %>% as.data.frame
gwas_coding$summary = summary$summary[match(gwas_coding$Symbol, summary$query)]

shinyList <- list(addi = addi , apoe2 = apoe2, apoe4 = apoe4)
save(shinyList, file = "Shiny/release_V3/data.rdt")

# intron

intron <- filter(vep.lod, grepl("intron_variant", Consequence))
(intron_gene <- intron$Symbol %>% unique)
intron_gene_gk <- myGK(intron_gene)

# non_coding_transcript_exon_variant

non_coding <- filter(vep.lod, grepl("non_coding_transcript_variant", Consequence))
non_coding <- filter(vep.lod, grepl("non_coding_transcript_exon_variant", Consequence))
non_coding <- non_coding[! duplicated(paste0(non_coding$UID, non_coding$Symbol)), ]
non_coding <- cbind(non_coding, gwas[match(non_coding$UID, gwas$UID), ])
non_coding <- non_coding[order(non_coding$LOD, decreasing = T), ]

non_coding <- non_coding[with(non_coding, abs(pSnp) > 0.2), ]

y = intron_gene_gk[[1]][1]
y[grepl("metabolic", y$Term), ]$Symbols %>% c %>% unique

# regulatory_region_variant

filter(vep.lod, grepl("regulatory_region_variant", Consequence))
ensembl_regulation = useMart(biomart="ENSEMBL_MART_FUNCGEN",host="www.ensembl.org",dataset="hsapiens_motif_feature")
