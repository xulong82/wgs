library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/Adsp")

load("./Manu/data.rdt")

mdata <- data$pheno
col.manual <- c("grey70", "firebrick1", "dodgerblue3", "gold1", "chartreuse3", "darkorchid2")
col.ad4 <- c("#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")

# --- PDF: Pie chart for AD and APOE

table <- table(mdata$AD2)  # AD status
pct <- round(table/sum(table)*100)
lbs <- names(table)
lbs <- paste(paste(lbs, pct), "%", sep = "")
pdf("./PDF/ad_pie.pdf", family = "Helvetica")
# png("./Shiny/www/adsp.png")
par(cex = 1.7, col = "grey30")
pie(table, labels = lbs, col = col.ad4, border = F)
dev.off()

table <- table(mdata$APOE)  # APOE
pct <- round(table/sum(table)*100)
lbs <- paste(names(table), ":", sep = "")
lbs <- paste(paste(lbs, pct), "%", sep = "")
pdf("./PDF/apoe_pie.pdf", family = "Helvetica")
par(cex = 1.7, col = "grey30")
pie(table, label = lbs, col = col.manual, border = F)
dev.off()

table <- table(mdata$Sex)  # Sex
pct <- round(table/sum(table)*100)
lbs <- c("Male", "Female")
lbs <- paste(paste(lbs, pct), "%", sep = "")
pdf("./PDF/sex_pie.pdf", family = "Helvetica")
par(cex = 1.7, col = "grey30")
pie(table, label = lbs, col = c("#1f78b4", "#fb9a99"), border = F)
dev.off()

# --- PDF: Barplot for samples per APOE

pdf("./Manu/apoe_bar.pdf", width = 4.5, height = 3, family = "Helvetica")

ggplot(mdata, aes(x = AD2, fill = APOE)) + geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("Genotype Frequency")  + coord_flip() +
  scale_fill_manual(values = col.manual,
                    name="APOE", breaks=c("33", "22", "23", "24", "34"), 
                    labels=c("3/3", "2/2", "2/3", "2/4", "3/4")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.key = element_blank()) 

dev.off()

# --- PDF: Histogram distribution for Age

pdf("./PDF/age_hist.pdf", width = 5, height = 7, family = "Helvetica")

ggplot(mdata, aes(x = Age)) + 
  geom_histogram(colour = "black", fill = "grey70", binwidth = 1) +
  theme_bw() + xlab("") + ylab("Case Number") +
  theme(panel.border = element_rect(size = 0.75, color = "grey30")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)) 

dev.off()

# --- PDF: Barplot for samples per Sex

pdf("./Manu/sex_bar.pdf", width = 4.5, height = 3, family = "Helvetica")

ggplot(mdata, aes(x = AD2, fill = Sex)) + geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("Sex Frequency") + coord_flip() +
  scale_fill_manual(values = c("#1f78b4", "#fb9a99"),
                    name="", breaks=c("0", "1"), labels=c("M", "F")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.key = element_blank()) 

dev.off()

# --- PDF: Boxplot for Age per AD status

pdf("./Manu/age_box.pdf", width = 4, height = 3, family = "Helvetica")

ggplot(mdata, aes(x = AD2, y = Age, fill = AD2)) + geom_boxplot(width = 0.75) +
# geom_hline(yintercept = 65, color = "darkred", size = 2, linetype = 2) + 
  scale_fill_manual(values = col.ad4) +
  theme_bw() + xlab("") + ylab("Age") + coord_flip() +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text = element_text(size = 12),
        legend.position = "none")

dev.off()

# --- PDF: Boxplot for Age by Sex per AD status

pdf("./PDF/age_sex_box.pdf", width = 4, height = 3, family = "Helvetica")
ggplot(mdata, aes(x = AD2, y = Age, fill = Sex)) + geom_boxplot(width = 0.75) +
  scale_fill_manual(values = c("#1f78b4", "#fb9a99"),
                    name="", breaks=c("0", "1"), labels=c("M", "F")) +
  theme_bw() + xlab("") + ylab("Age") + # coord_flip() +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12))
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
