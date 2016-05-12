library(rstan)
library(reshape)
library(dplyr)
library(xlsx)
library(stargazer)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")
load("./Manu/data.rdt") 
load("./data/kinship.rdt")
load("./Manu/model.rdt")

u = model$optimizing$par
u = u[paste0("u[", 1:576, "]")]

hist(u)

pheno = data$pheno
pheno$random = u

options(stringsAsFactors = F)
meta_extra <- read.xlsx("./docs/adsp_subject_and_experiment_metadata.xlsx", startRow = 2, sheetIndex = 1)

pheno$Ethnicity <- meta_extra$Ethnicity[match(pheno$SRR, meta_extra$SRR)]
pheno$Ethnicity[pheno$Ethnicity == "0"] = "Not-Latino"
pheno$Ethnicity[pheno$Ethnicity == "1"] = "Latino"
pheno$Ethnicity[pheno$Ethnicity == "null"] = "Null"

pheno$Race[pheno$Race == "4"] = "Black"
pheno$Race[pheno$Race == "5"] = "White"
pheno$Race[pheno$Race == "6"] = "Other"
pheno$Race[pheno$Race == "null"] = "Null"

pheno$Race = factor(pheno$Race, levels = c("Black", "White", "Other", "Null"))

table(pheno$Race, pheno$Ethnicity)
stargazer(table(pheno$Race, pheno$Ethnicity) %>% as.data.frame.matrix, summary = F)

mycol <- c("grey10", "dodgerblue3", "firebrick1", "chartreuse3")

pdf("./Manu/random.pdf", width = 7, height = 12, family = "Helvetica")

ggplot(pheno, aes(x = random, y = Family.ID)) + 
  geom_line() + geom_point(aes(color = Race, shape = Ethnicity)) +
  theme_bw() + xlab("Random Effect") + ylab("Family ID") +
  scale_color_manual(values = mycol) +
  theme(panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.text = element_text(size = 7),
        axis.title= element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key = element_blank()) 

dev.off()

dt = data.frame(u = u)

pdf("./Pdf/random.pdf", width = 6, height = 5, family = "Helvetica")

ggplot(dt, aes(x = u)) +
  geom_histogram(fill = "dodgerblue3", color = "white", size = .5) +
  theme_bw() + xlab("") + ylab("Count") +
  theme(panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.text = element_text(size = 15),
        axis.title= element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_blank(), 
        legend.key = element_blank()) 

dev.off()

k = kinship$autosome

u2 = matrix(0, nrow = 576, ncol = 576)

cor(c(k), c(u2))

for (i in 2:576)
  for (j in 1:(i-1)) 
    u2[i, j] = u2[j, i] = abs(u[i] - u[j])
