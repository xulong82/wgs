chr <- commandArgs(TRUE)  # 1-22, X, Y

# geno <- list()
# gInf <- NULL
# load("~/Dropbox/ADSP/R/mdata.rdt")
# map <- data.frame(row.names = c("0/0", "0/1", "1/1", "1/2"), g = c(0, 1, 2, 3))
# name <- list.files(path = "/data/xwang/ADSP/VCF", pattern = "*.vcf")
# 
# # MAKE GENOTYPE MATRIX & DEFINE VARIANTS
# for (i in 1:length(name)) {
#   filepath <- file.path("/data/xwang/ADSP/VCF", name[i])
#   geno0 <- read.delim(filepath, sep = "", stringsAsFactors = F)
# 
#   geno0 <- geno0[geno0$CHR == chr, ]
#   geno0$UID <- paste(geno0$CHR, geno0$POS, sep = "-")
#   geno0$FID <- apply(geno0[c(7, 4:5)], 1, function (x) paste(x, collapse = "-"))
#   geno0$GTN <- map[geno0$GT, "g"]
#   rownames(geno0) <- NULL
#   geno[[i]] <- geno0[c("UID", "FID", "GTN")]
# 
#   gInf0 <- geno0[! geno0$FID %in% gInf$FID, ]
#   gInf <- rbind(gInf, gInf0[c("UID", "CHR", "POS", "ID", "REF", "ALT", "FID")])
# }
# 
# names(geno) <- gsub(".vcf", "", name)
# geno <- geno[mdata$SRR]  
# rownames(gInf) <- NULL
#   
# setwd("/data/xwang/ADSP/geno")
# fname1 <- paste("gInf_chr", chr, ".rdt", sep = "")
# fname2 <- paste("geno_chr", chr, ".rdt", sep = "")
# 
# save(gInf, file = fname1)
# save(geno, file = fname2)
#   
# load(fname1)
# load(fname2)

# CLASSIFY VARIANTS
# base4 <- c("A", "T", "C", "G")
# myidx <- gInf$REF %in% base4 & gInf$ALT %in% base4
# snpInf <- gInf[myidx, ]
# mixInf <- gInf[! myidx, ]
# polyInf <- mixInf[grepl("," , mixInf$ALT), ]
# 
# snpInf <- snpInf[! snpInf$UID %in% mixInf$UID, ]
# pSnpId <- snpInf$UID[duplicated(snpInf$UID)]
# snpInf <- snpInf[! snpInf$UID %in% pSnpId, ]  # di-SNP
# 
# myidx <- polyInf$REF %in% base4 & nchar(polyInf$ALT) == 3
# pSnpId <- unique(c(pSnpId, polyInf$UID[myidx]))
# pSnpInf <- gInf[gInf$UID %in% pSnpId, ]  # poly-SNP
# 
# indInf <- mixInf[! mixInf$UID %in% polyInf$UID, ]
# pIndId <- indInf$UID[duplicated(indInf$UID)]
# indInf <- indInf[! indInf$UID %in% pIndId, ]  # Di-INDEL
# 
# pIndId <- unique(c(pIndId, polyInf$UID[! nchar(polyInf$ALT) == 3]))
# pIndInf <- gInf[gInf$UID %in% pIndId, ]  # Poly-INDEL

# length(unique(gInf$UID)) 
# length(intersect(pSnpId, pIndId))  # loci with both polyINDEL and polySNP characteristic
# nrow(snpInf) + length(unique(pSnpInf$UID)) + nrow(indInf) + length(unique(pIndInf$UID))

# gMix <- list()
# gMix$snpInf <- snpInf
# gMix$pSnpInf <- pSnpInf
# gMix$indInf <- indInf
# gMix$pIndInf <- pIndInf
  
# setwd("/data/xwang/ADSP/gMix")
# fname1 <- paste("gMix_chr", chr, ".rdt", sep = "")

# save(gMix, file = fname1)

setwd("/data/xwang/ADSP")
load(paste0("geno/geno_chr", chr, ".rdt"))  # geno
load(paste0("gMix/gMix_chr", chr, ".rdt"))  # gMix

n.sample <- length(geno)
id.sample <- names(geno)

# bi-Var-only genotype
biInf <- rbind(gMix$snpInf, gMix$indInf)
biInf <- biInf[order(biInf$POS), ]
n.var <- nrow(biInf)
id.var <- biInf$UID

biGeno <- matrix(0, n.var, n.sample, dimnames = list(id.var, id.sample))
for (i in 1:n.sample) {
  geno0 <- geno[[i]]
  geno0 <- geno0[geno0$UID %in% id.var, ]
  biGeno[match(geno0$UID, id.var), i] <- geno0[, "GTN"]
}

save(biInf, file = paste0("biVar/biInf_chr", chr, ".rdt"))
save(biGeno, file = paste0("biVar/biGeno_chr", chr, ".rdt"))
write.table(biInf, file = paste0("biVar/biInf_chr", chr, ".txt"), quote = F, row.names = F, col.names = F)
write.table(biGeno, file = paste0("biVar/biGeno_chr", chr, ".txt"), row.names = F, col.names = F)

chr <- paste0("chr", c(1:22, "X", "Y"))
varAll_sample <- sapply(chr, function(x) { cat(x, "\n")
  load(paste0("geno/geno_", x, ".rdt"))
  sapply(geno, nrow)
})  # SAVED

