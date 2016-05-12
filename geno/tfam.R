# Copyright: Xulong Wang (xulong.wang@jax.org)

# ADSP.TFAM FILE ---
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype

load("~/Dropbox/GitHub/Adsp/data/mdata.rdt")
adsp.tfam <- data.frame(row.names = 1:nrow(mdata))
adsp.tfam$Family.ID <- as.numeric(as.factor(mdata$Family.ID))
adsp.tfam$Individual.ID <- 1:nrow(mdata)
adsp.tfam$Paternal.ID <- rep(0, nrow(mdata))
adsp.tfam$Maternal.ID <- rep(0, nrow(mdata))
adsp.tfam$Sex <- as.numeric(mdata$Sex)
adsp.tfam$Phenotype <- mdata$AD1

setwd("/data/xwang/ADSP/Plink")
write.table(adsp.tfam, file = "adsp.tfam", row.names = F, col.names = F)

adsp.cov <- data.frame(row.names = 1:nrow(mdata))
adsp.cov$Family.ID <- as.numeric(as.factor(mdata$Family.ID))
adsp.cov$Individual.ID <- 1:nrow(mdata)
adsp.cov$Age <- mdata$Age
adsp.cov$Sex <- as.numeric(mdata$Sex) - 1

setwd("/data/xwang/Adsp/Plink")
write.table(adsp.cov, file = "adsp_cov.txt", row.names = F, col.names = F)

