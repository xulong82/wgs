library(scales)
library(ggplot2)
library(dplyr)

rm(list = ls())
setwd("~/Dropbox/GitHub/wgs")

load("./new/mdata.rdt")
mdata$APOE = factor(mdata$APOE, levels = c("33", "22", "23", "24", "34"))

col4 <- c("#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")
col6 <- c("grey70", "firebrick1", "dodgerblue3", "gold1", "chartreuse3", "darkorchid2")

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    plot.title=element_text(size=14, face="bold")
)

# --- PDF: Pie chart for AD and APOE

df <- as.data.frame(table(mdata$AD1))
names(df) = c("name", "tbl")
df = mutate(df, perc = tbl / sum(tbl), label_pos = cumsum(perc) - perc / 2,  perc_text = paste0(round(perc * 100), "%"))

pdf("./Manu/AD.pdf", width = 4, height = 4)

ggplot(df, aes(x = "", y = perc, fill = name)) + 
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) +
  blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = col4) +
  geom_text(aes(y = label_pos, label = perc_text), size = 3)

dev.off()

tbl <- table(mdata$AD1) # AD
pct <- round(tbl/sum(tbl)*100)
lbs <- paste(names(tbl), ":", sep = "")
lbs <- paste(paste(lbs, pct), "%", sep = "")

pdf("./Manu/AD.pdf")
par(mar = c(5, 4, 4, 4), cex = 1.8, col = "grey30")
pie(tbl, labels = lbs, col = col4, border = F, init.angle = 45)
dev.off()

tbl <- table(mdata$APOE) # APOE
names(tbl) = c("e3/e3", "e2/e2", "e2/e3", "e2/e4", "e3/e4")
pct <- round(tbl/sum(tbl)*100)
lbs <- paste(names(tbl), ":", sep = "")
lbs <- paste(paste(lbs, pct), "%", sep = "")

pdf("./Manu/APOE.pdf")
par(cex = 1.8, col = "grey30")
pie(tbl, label = lbs, col = col6, border = F)
dev.off()

tbl <- table(mdata$Sex)  # Sex
pct <- round(tbl/sum(tbl)*100)
lbs <- c("Male", "Female")
lbs <- paste(paste(lbs, pct), "%", sep = "")

pdf("./Manu/Sex.pdf")
par(cex = 1.8, col = "grey30")
pie(tbl, label = lbs, col = c("#1f78b4", "#fb9a99"), border = F)
dev.off()

# --- PDF: Barplot for samples per APOE

pdf("./Manu/APOE_bar.pdf", width = 4, height = 2)

ggplot(mdata, aes(x = AD1, fill = APOE)) + geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("Genotype Frequency") + coord_flip() +
  scale_fill_manual(values = col6, name="APOE", breaks=c("33", "22", "23", "24", "34"), 
                    labels=c("e3/e3", "e2/e2", "e2/e3", "e2/e4", "e3/e4")) +
  theme(legend.title = element_blank(), legend.key = element_blank()) 

dev.off()

# --- PDF: Barplot for samples per Sex

pdf("./Manu/Sex_bar.pdf", width = 4, height = 2)

ggplot(mdata, aes(x = AD1, fill = as.factor(Sex))) + geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("Sex Frequency") + coord_flip() +
  scale_fill_manual(values = c("#1f78b4", "#fb9a99"), name="", breaks=c("0", "1"), labels=c("M", "F")) +
  theme(legend.title = element_blank(), legend.key = element_blank()) 

dev.off()

# --- PDF: Boxplot for Age per AD status

pdf("./Manu/Age.pdf", width = 3, height = 2)

ggplot(mdata, aes(x = AD1, y = Age, fill = AD1)) + geom_boxplot(width = 0.65) +
  scale_fill_manual(values = col4) + theme_bw() + xlab("") + ylab("Age") + coord_flip() +
  theme(legend.position = "none")

dev.off()

# --- PDF: Barplot for samples per family

pdf("./PDF/family_bar.pdf", width = 12, height = 8, family = "Helvetica")
ggplot(mdata, aes(x = Family.ID, fill = AD2)) + geom_bar(width = 0.65) +
  theme_bw() + xlab("") + ylab("Case Number") +
  scale_fill_manual(values = col.ad4) +
  theme(panel.border = element_rect(size = 0.75, color = "grey30")) +
  theme(axis.text.x = element_text(size = 6, angle = -90),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15, vjust = 1)) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()
