library(lme4)
library(ordinal)
library(lmtest)
library(parallel)

load("/data/xwang/ADSP/R/mdata1.rdt")
load("/data/xwang/ADSP/R/vdata.rdt")

cat(date(), "--- create data for parallel computing \n")
geno <- geno[sample(1:nrow(geno), 1e6), ]  # random 1M snps

n.core <- 20  # maxLod in 1e6/20 permutated snps
n.snp <- nrow(geno)
n.sp <- ncol(geno)
unit <- round(n.snp / n.core) 
start <- ((1:n.core) - 1) * unit + 1
end <- (1:n.core) * unit
end[n.core] <- n.snp

genoList <- list()
for (i in 1:n.core) 
  genoList[[i]] <- geno[start[i]:end[i], ]

cat(date(), "--- multilevel linear regression \n")
data1 <- mdata[, c("AD1", "Sex", "Age", "APOE", "Family.ID")]
data1$AD2 = factor(data1$AD1, levels = c("0", "0.25", "0.5", "1"), ordered = T)

# model1 <- lm(formula = AD1 ~ Sex + Age, data1)
# model2 <- lm(formula = AD1 ~ Sex + Age + APOE, data1)
# model3 <- lmer(formula = AD1 ~ Sex + Age + (1 | Family.ID), data1, REML = FALSE)
# model4 <- lmer(formula = AD1 ~ Sex + Age + APOE + (1 | Family.ID), data1, REML = FALSE)
 
model1 <- clm(formula = AD2 ~ Sex + Age, data = data1)
model2 <- clm(formula = AD2 ~ Sex + Age + APOE, data = data1)
model3 <- clmm(formula = AD2 ~ Sex + Age + (1 | Family.ID), data1)
model4 <- clmm(formula = AD2 ~ Sex + Age + APOE + (1 | Family.ID), data1)

loglik10 <- logLik(model1)
loglik20 <- logLik(model2)
loglik30 <- logLik(model3)
loglik40 <- logLik(model4)
loglik0 <- c(loglik10, loglik20, loglik30, loglik40)

mixLinear <- function (x) {
  n1 <- nrow(x)
  maxLogLik = rep(-1e3, 4)
  for (i in 1:n1) {
    data2 <- cbind(data1, SNP = sample(x[i, ]))

#   model1 <- lm(formula = AD1 ~ Sex + Age + SNP, data2)
#   model2 <- lm(formula = AD1 ~ Sex + Age + APOE + SNP, data2)
#   model3 <- lmer(formula = AD1 ~ Sex + Age + SNP + (1 | Family.ID), data2, REML = FALSE)
#   model4 <- lmer(formula = AD1 ~ Sex + Age + APOE + SNP + (1 | Family.ID), data2, REML = FALSE)

    model1 <- clm(formula = AD2 ~ Sex + Age + SNP, data = data2)
    model2 <- clm(formula = AD2 ~ Sex + Age + APOE + SNP, data = data2)
    model3 <- clmm(formula = AD2 ~ Sex + Age + SNP + (1 | Family.ID), data2)
    model4 <- clmm(formula = AD2 ~ Sex + Age + APOE + SNP + (1 | Family.ID), data2)

    if (logLik(model1) > maxLogLik[1]) maxLogLik[1] <- logLik(model1)
    if (logLik(model2) > maxLogLik[2]) maxLogLik[2] <- logLik(model2)
    if (logLik(model3) > maxLogLik[3]) maxLogLik[3] <- logLik(model3)
    if (logLik(model4) > maxLogLik[4]) maxLogLik[4] <- logLik(model4)
  }
  return(maxLogLik)
}

cat(date(), "--- multiple cores run \n")
maxLogLik <- matrix(nrow = n.core * 4, ncol = 4)
for (i in 1:4) {  # 20 by 4 permutations
  y <- mclapply(genoList, mixLinear, mc.cores = 20)
  start <- (i - 1) * 20 + 1
  end <- i * 20
  maxLogLik[start:end, ] <- matrix(unlist(y), ncol = 4, byrow = T)
}
cat(date(), "--- fi \n")

arg <- commandArgs(TRUE)
fname <- paste("maxLogLik", arg, sep = "")
assign(fname, sweep(maxLogLik, 2, loglik0, "-"))
# save(list = fname, file = paste(paste("/data/xwang/ADSP/R/permut1", fname, sep = "/"), "rdt", sep = "."))
save(list = fname, file = paste(paste("/data/xwang/ADSP/R/permut2", fname, sep = "/"), "rdt", sep = "."))

library(GGally)
library(ggplot2)

load("/data/xwang/ADSP/dbSNP/dbSNP137.rdt")

cat("--- LM logLik data \n")
load("/data/xwang/ADSP/R/logLik_lm.rdt")

cat("--- merge logLik data \n")
fname1 <- paste("logLikClm", 1:10, sep = "")
fname2 <- paste(fname1, "rdt", sep = ".")
load(paste("/data/xwang/ADSP/R/Clm_novo", fname2[1], sep = "/"))
logLik <- get(fname1[1])
for (i in 2:10) {
  load(paste("/data/xwang/ADSP/R/Clm_novo", fname2[i], sep = "/"))
  logLik <- rbind(logLik, get(fname1[i]))
}
sscan <- as.data.frame(logLik)
sscan$CHROM <- gsub("-.*", "", rownames(logLik))
sscan$POS <- gsub("^.*-", "", rownames(logLik))
sscan$POS <- as.numeric(sscan$POS)
sscan$ID <- rownames(logLik)
save(sscan, file = "/data/xwang/ADSP/R/sscan_clm_novo.rdt")

table(rownames(logLik) %in% dbSNP137$ID)
#--- caveat: all TRUE or problem
snp.inf <- dbSNP137[match(rownames(logLik), dbSNP137$ID), ]
sscan <- cbind(snp.inf, logLik)
save(sscan, file = "/data/xwang/ADSP/R/sscan_lm_rs.rdt")
save(sscan, file = "/data/xwang/ADSP/R/sscan_clm_rs.rdt")

sscan1 <- sscan[apply(sscan[, 4:7], 1, function (x) {max(x) > 3}), ]
save(sscan1, file = "/data/xwang/ADSP/R/sscan_lm1.rdt")

sscan1 <- sscan[apply(sscan[, 4:7], 1, function (x) {max(x) > 2}), ]
save(sscan1, file = "/data/xwang/ADSP/R/sscan_clm1.rdt")

png(file = "/data/xwang/ADSP/R/ggpair_lm.png", width = 480, height = 480)
png(file = "/data/xwang/ADSP/R/ggpair_clm.png", width = 480, height = 480)
ggpairs(sscan1, columns = 4:7,
 	upper= list(continuous = "cor", params = c(corSize = 10)), 
        lower = list(continuous = "smooth", params = c(color = "red")),
	diag = list(continuous = "density"),
	axisLabels = "show")
dev.off()

# --- merge the maxLogLik data \n")
fname1 <- paste("maxLogLik", 1:20, sep = "")
fname2 <- paste(fname1, "rdt", sep = ".")
fname3 <- paste("/data/xwang/ADSP/R/permut2", fname2, sep = "/")
load(fname3[1])
maxLogLik <- get(fname1[1])
for (i in 2:20) {
  if (file.exists(fname3[i])) {
    load(fname3[i])
    maxLogLik <- rbind(maxLogLik, get(fname1[i]))
  }
}
save(maxLogLik, file = "~/Dropbox/ADSP/R/permut_lm.rdt")
save(maxLogLik, file = "~/Dropbox/ADSP/R/permut_clm.rdt")

rm(list = ls())

load("~/Dropbox/ADSP/R/hitRs_clm.rdt")
snpIdpool <- hit.rs$m2

snpId <- "rs34572242"  # ADSP Locus

load("/data/xwang/ADSP/R/vdata_rs.rdt")
geno <- geno[snpIdpool, ]

load("/data/xwang/ADSP/R/mdata1.rdt")
data <- cbind(mdata, t(geno))

for (snpId in snpIdpool) {
  cat(snpId, "\n")
# data$geno <- data[, "rs6689933"]
  data$geno <- data[, snpId]
  
  pctg0 <- rep(0, 4)
  names(pctg0) <- c("0", "0.25", "0.5", "1")
  pctg1 <- pctg2 <- pctg0

  pctg <- table(data[data$geno == 0, ]$AD1) / sum(data$geno == 0)
  pctg0[names(pctg)] <- pctg

  pctg <- table(data[data$geno == 1, ]$AD1) / sum(data$geno == 1)
  pctg1[names(pctg)] <- pctg

  pctg <- table(data[data$geno == 2, ]$AD1) / sum(data$geno == 2)
  pctg2[names(pctg)] <- pctg

  graph.dt <- data.frame(value = c(pctg0, pctg1, pctg2), 
                         AD = rep(c("No", "Possible", "Probable", "Definite"), 3), 
                         geno = rep(c("0", "1", "2"), each = 4))
  graph.dt$AD <- factor(graph.dt$AD, levels = c("No", "Possible", "Probable", "Definite"))
  
  filepath <- paste(paste("~/Dropbox/ADSP/group", snpId, sep = "/"), "pdf", sep = ".")
  pdf(filepath, width = 5)
  ggplot(graph.dt) + geom_line(aes(x = as.numeric(geno), y = value, colour = AD), size = 2) +
    scale_x_continuous(breaks = c(1:3), labels = c("0/0", "0/1", "1/1")) +
    theme_bw() + xlab("") + ylab("Percentage %") +
    theme(panel.border = element_rect(size = 1, color = "black")) +
    theme(axis.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 15, face = "bold")) +
    theme(legend.position = "top", legend.direction = "horizontal", 
          legend.text = element_text(size = 12, face = "bold"),
          legend.title = element_blank(), legend.key = element_blank()) 
  dev.off()
}

# ggplot(data, aes(as.factor(geno), fill = AD2)) + geom_bar()
# ggplot(data, aes(as.factor(geno), fill = AD2)) + geom_bar(position = "fill")
# ggplot(data, aes(AD2, fill = as.factor(geno))) + geom_bar()
# ggplot(data, aes(AD2, fill = as.factor(geno))) + geom_bar(position = "fill")

# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: GWAS logLik graph

library(ggplot2)
library(GGally)
library(xtable)

rm(list = ls())

#--- LogLik in single scan
load("/data/xwang/ADSP/R/sscan_clm_rs.rdt")  # rs-prefixed
sscan1 <- sscan

load("/data/xwang/ADSP/R/sscan_clm_novo.rdt")  # de-novo
sscan$CHROM <- paste("chr", sscan$CHROM, sep = "")
sscan <- rbind(sscan1, sscan[, c("CHROM", "POS", "ID", "m1", "m2", "m3", "m4")])
rm(sscan1)

source("~/Dropbox/ADSP/R/peaks.R")  # automatic peak detection

chr <- read.delim("~/Dropbox/X/chromInfo.hg19.txt", header = F, stringsAsFactors = F)
chrlen <- cumsum(as.numeric(chr$V2)) * 1e-6  # mbp unit
names(chrlen) <- c(paste("chr", c(1:22, "X", "Y"), sep = ""))
chrmid <- diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])

#--- genomic logLik graph
sscan$CHROM1 <- sscan$CHROM
sscan$CHROM1 <- gsub("chrX", "chr23", sscan$CHROM1)
sscan$CHROM1 <- gsub("chrY", "chr24", sscan$CHROM1)
sscan$CHROM1 <- as.numeric(gsub("chr", "", sscan$CHROM1))
sscan$POS1 <- c(0, chrlen)[sscan$CHROM1] + sscan$POS * 1e-6 

sscan$CHROM <- factor(sscan$CHROM, levels = c(paste("chr", c(1:22, "X", "Y"), sep = "")))
sscan$POS <- sscan$POS * 1e-6  # mb

#--- manhatton
sscan$odds <- rep(0, nrow(sscan))
sscan$odds[sscan$CHROM1 %% 2 == 1] <- 1
png("~/Dropbox/ADSP/logLik/manhattan_clm_m4.png", width = 2048, height = 1024)
  ggplot(sscan) + geom_point(aes(x = POS1, y = m4, color = as.factor(odds)), size = 2) +
#   geom_rect(data = data.frame(x1 = chrlen[2*(1:12) - 1], x2 = chrlen[2*(1:12)]), 
#             aes(xmin = x1, xmax = x2, ymin = -Inf, ymax = +Inf),
#             alpha = 0.2, fill = "lightblue") +
  theme_bw() + xlab("") + ylab("Log Likelihood Ratio") +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) + 
  scale_x_continuous(breaks = chrmid, labels = names(chrlen)) + ylim(0, 18) +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text.x = element_text(size = 25, angle = -90, face = "bold", 
                      colour = rep(c("#377eb8", "#e41a1c"))),
        axis.text.y = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 30, face = "bold")) +
  theme(legend.position = "none")
dev.off()

#--- chromsome by chromsome
# sscan1 <- sscan[apply(sscan[, 4:7], 1, function (x) {max(x) > 6}), ]
load("~/Dropbox/ADSP/R/sscan_clm1.rdt")  # test purpose
sscan1$POS <- sscan1$POS * 1e-6
sscan1$hit.m1 <- rep("no", nrow(sscan1))
sscan1$hit.m4 <- sscan1$hit.m3 <- sscan1$hit.m2 <- sscan1$hit.m1
load("~/Dropbox/ADSP/R/hitRs_clm.rdt")
sscan1$hit.m1[sscan1$ID %in% hit.rs[["m1"]]] = "yes"
sscan1$hit.m2[sscan1$ID %in% hit.rs[["m2"]]] = "yes"
sscan1$hit.m3[sscan1$ID %in% hit.rs[["m3"]]] = "yes"
sscan1$hit.m4[sscan1$ID %in% hit.rs[["m4"]]] = "yes"

pdf(file = "~/Dropbox/ADSP/logLik/logLik_clm1.pdf", width = 15, height = 10)
  ggplot(sscan1, aes(x = POS, y = m1)) + geom_point(aes(color = hit.m1, size = hit.m1)) + 
pdf(file = "~/Dropbox/ADSP/logLik/logLik_clm2.pdf", width = 15, height = 10)
  ggplot(sscan1, aes(x = POS, y = m2)) + geom_point(aes(color = hit.m2, size = hit.m2)) + 
pdf(file = "~/Dropbox/ADSP/logLik/logLik_clm3.pdf", width = 15, height = 10)
  ggplot(sscan1, aes(x = POS, y = m3)) + geom_point(aes(color = hit.m3, size = hit.m3)) + 
pdf(file = "~/Dropbox/ADSP/logLik/logLik_clm4.pdf", width = 15, height = 10)
  ggplot(sscan1, aes(x = POS, y = m4)) + geom_point(aes(color = hit.m4, size = hit.m4)) + 
  #----------------------------------------------------------------------------------------
  scale_size_discrete(range = c(0.5, 2)) + scale_color_manual(values = c("lightblue", "red")) +
  # geom_hline(yintercept = 12, colour = "gold", size = .1) +
  facet_grid(CHROM ~ .) + theme_bw() + xlab("Genome (Mbp)") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, face = "bold", color = "grey30"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_blank()) + theme(legend.position = "none")
dev.off()

table.tex <- sscan1[sscan1$ID %in% hit.rs[["m4"]], ]
table.tex <- table.tex[, c("ID", "CHROM", "POS", "m1", "m2", "m3", "m4")]
rownames(table.tex) <- c(1:nrow(table.tex))
table.tex$POS <- as.integer(table.tex$POS)
print(xtable(table.tex, caption = "model4"))
# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: cutoffs with permutations

rm(list = ls())

#--- maxLogLik distribution in permutations
load("~/Dropbox/ADSP/R/permut_lm.rdt")  # linear models
load("~/Dropbox/ADSP/R/permut_clm.rdt")  # cumulative link models
colnames(maxLogLik) <- c("model1", "model2", "model3", "model4")
maxLogLik.df <- as.data.frame(maxLogLik)

maxLogLik.df$hit.m1 <- rep("no", nrow(maxLogLik.df))
maxLogLik.df$hit.m4 <- maxLogLik.df$hit.m3 <- maxLogLik.df$hit.m2 <- maxLogLik.df$hit.m1

cut.model1 <- quantile(maxLogLik.df$model1, 0.95)
cut.model2 <- quantile(maxLogLik.df$model2, 0.95)
cut.model3 <- quantile(maxLogLik.df$model3, 0.95)
cut.model4 <- quantile(maxLogLik.df$model4, 0.95)

maxLogLik.df$hit.m1[maxLogLik.df$model1 > cut.model1] = "yes"
maxLogLik.df$hit.m2[maxLogLik.df$model2 > cut.model2] = "yes"
maxLogLik.df$hit.m3[maxLogLik.df$model3 > cut.model3] = "yes"
maxLogLik.df$hit.m4[maxLogLik.df$model4 > cut.model4] = "yes"

permut.dt <- data.frame(value = c(as.matrix(maxLogLik.df[, 1:4])), 
                        yes = factor(c(as.matrix(maxLogLik.df[, 5:8])), levels = c("yes", "no")),
                        model = rep(c("model1", "model2", "model3", "model4"), each = nrow(maxLogLik)))

pdf("~/Dropbox/ADSP/Models/permut_lm.pdf")
ggplot(permut.dt, aes(x = value)) + geom_histogram(aes(fill = yes), binwidth = .1) +
  facet_grid(model ~ .) +
  scale_fill_manual(values = c("#1f78b4", "#fb9a99")) +
  theme_bw() + xlab("") + ylab("Max LRT Count") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

