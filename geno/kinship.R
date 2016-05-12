rm(list = ls())

setwd("/data/xwang/ADSP/Kinship")
load("~/Dropbox/GitHub/Adsp/data/mdata.rdt")

kin <- read.delim("autosome.kin")
kin0 <- read.delim("autosome.kin0")

idx <- rbind(kin[, c("ID1", "ID2")], kin0[, c("ID1", "ID2")])
idx <- t(apply(idx, 1, function(x) sort(x)))

kin23 <- c(kin$Kinship, kin0$Kinship)

KS0 <- matrix(0, nrow(mdata), nrow(mdata))
colnames(KS0) <- rownames(KS0) <- mdata$SRR

KS <- KS0
for (j in 1:nrow(idx)) KS[idx[j, 2], idx[j, 1]] <- kin23[j]
diag(KS) <- 1
KS[upper.tri(KS)] <- t(KS)[upper.tri(KS)]
kinship <- list(autosome = KS)

for (x in paste0("no_chr", 1:22)) { cat(x, "\n")
  kin <- read.delim(paste(x, "kin", sep = "."))
  kin0 <- read.delim(paste(x, "kin0", sep = "."))
  kin230 <- c(kin$Kinship, kin0$Kinship)
  kin23 <- cbind(kin23, kin230)
   
  KS <- KS0
  for (j in 1:nrow(idx)) KS[idx[j, 2], idx[j, 1]] <- kin230[j]
  diag(KS) <- 1
  KS[upper.tri(KS)] <- t(KS)[upper.tri(KS)]
  kinship[[x]] <- KS
}

colnames(kin23) <- c("autosome", paste0("no_chr", 1:22))
kinship$kin23 <- kin23

# --- LOCAL ---
setwd("~/Dropbox/GitHub/Adsp")
save(kinship, file = "./data/kinship.rdt")
