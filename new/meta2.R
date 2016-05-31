rm(list = ls())

load("~/Dropbox/GitHub/wgs/new/mdata.rdt")

SRR <- mdata$SRR
Family <- mdata$Family.ID

meta = lapply(1:22, function(chr) {

cat(chr, "\n")

load(paste0("/data/xwang/adsp3/R/chr", chr, ".rdt"))
geno = geno[, mdata$ADSP.Sample.ID]

meta = data.frame(ID = rownames(geno))
meta$CHR = gsub(":.*", "", meta$ID)
meta$POS = gsub("^.*:(.*)_.*", "\\1", meta$ID)
meta$REF = gsub("^.*_(.*)/.*", "\\1", meta$ID)
meta$ALT = gsub(".*/", "", meta$ID)

meta$MAF = rowSums(geno) / ncol(geno) / 2

meta$N.SRR <- apply(geno, 1, function (x) sum(x > 0))
meta$N.Family <- apply(geno, 1, function (x) length(unique(Family[x > 0])))

meta

})

meta = do.call(rbind, meta)

save(meta, file = "/data/xwang/adsp3/meta.rdt")

