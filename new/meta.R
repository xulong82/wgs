#------- Note: make tfam file

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")

options(stringsAsFactors = F)
mdata <- read.delim("./docs/adsp_subject_and_experiment_metadata.txt")
table(duplicated(mdata$ADSP.Sample.ID))

mdata = mdata[mdata$AD %in% c(0, 1, 2, 3), ]
mdata = mdata[mdata$Age != "null", ]
mdata$Age = as.numeric(mdata$Age)

mdata$AD1 = "No"
mdata$AD1[mdata$AD == 3] = "Possible"
mdata$AD1[mdata$AD == 2] = "Probable"
mdata$AD1[mdata$AD == 1] = "Definite"

mdata$AD1 = factor(mdata$AD1, levels = c("No", "Possible", "Probable", "Definite"))

tfam <- read.table("/data/xwang/adsp3/plink/wgs.tfam")
tfam <- tfam[tfam$V1 %in% mdata$ADSP.Sample.ID, ]

mdata <- mdata[match(tfam$V1, mdata$ADSP.Sample.ID), ] # 570 only
all(tfam$V1 == mdata$ADSP.Sample.ID)

save(mdata, file = "./new/mdata.rdt")

# Plink TFAM FILE ---
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype (0=missing; 1=unaffected; 2=affected)

tfam$V5 = mdata$Sex + 1
tfam$V6 <- as.numeric(mdata$AD1) 

write.table(tfam, file = "/data/xwang/adsp3/plink/wgs2.tfam", row.names = F, col.names = F, quote = F)

