rm(list = ls())

# MAKE VCF 4 VEP

setwd("/data/xwang/ADSP")
chr <- c(1:22, "X", "Y")

biVcf <- lapply(chr, function(x) { cat(x)
  load(paste0("Meta/meta_chr", x, ".rdt"))
  meta = meta[meta$MAF > 0.01 & meta$HET < 0.99, ]
  out = meta[c("CHR", "POS", "UID", "REF", "ALT")]
  out$QUAL = out$FILTER = out$INFO = "."
  colnames(out)[1] = "#CHR"
  write.table(out, file = paste0("VEP/chr", x, ".vcf"), quote = F, row.names = F)
})

# PARSE VEP output 

rm(list = ls())
setwd("/data/xwang/ADSP/VEP")
chr <- paste("chr", c(1:22, "X", "Y"), sep = "")

vepList <- lapply(chr, function(x) { cat(x, "\n")
  vep1 <- read.delim(paste0(x, ".txt"), stringsAsFactors = F, comment.char = ".")
  vep1$Symbol <- vep1$Distance <- vep1$Biotype <- "-"

  vep1$Symbol[grep("SYMBOL", vep1$Extra)] <- 
    gsub("^.*SYMBOL=(.*);SYMBOL.*", "\\1", vep1$Extra[grep("SYMBOL", vep1$Extra)])
  vep1$Distance[grep("DISTANCE", vep1$Extra)] <- 
    gsub("^.*DISTANCE=(.*);STRAND.*", "\\1", vep1$Extra[grep("DISTANCE", vep1$Extra)])
  vep1$Biotype[grep("BIOTYPE", vep1$Extra)] <- 
    gsub("^.*BIOTYPE=(.*)", "\\1", vep1$Extra[grep("BIOTYPE", vep1$Extra)])
  vep1$Biotype <- gsub(";.*", "", vep1$Biotype)
  vep1$Extra <- NULL
  vep1
})

names(vepList) <- chr
save(vepList, file = "vepList.rdt")

vepStat <- lapply(vepList, function(x) {
  table(x$Consequence)
}) # SAVED

vep_UID_exon <- lapply(vepList, function(x) {
  x$X.Uploaded_variation[grep("exon", x$Consequence)] %>% unique
}) # NOT SAVED

vep_UID_intron <- lapply(vepList, function(x) {
  x$X.Uploaded_variation[grep("intron", x$Consequence)] %>% unique
}) # NOT SAVED

