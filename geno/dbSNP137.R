# --- Create R data for dbSNP137 database

dbSNP137 <- read.delim("/data/xwang/ADSP/dbSNP/dbSNP137Hg19.vcf", stringsAsFactors = F)
save(dbSNP137, file = "/data/xwang/ADSP/dbSNP/dbSNP137.rdt")

