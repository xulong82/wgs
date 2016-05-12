library(xlsx)
library(ggvis)

rm(list = ls())
hpc <- "/data/xwang/ADSP"
github <- "~/Dropbox/GitHub/Adsp"

setwd(github)
load("./data/hg19.rdt")
hg19 <- hg19[, c("NAME", "CHR", "STRAND", "START", "END", "POS")]
hg19$START <- hg19$START * 1e-3
hg19$END <- hg19$END * 1e-3
hg19$POS <- hg19$POS * 1e-3

chr <- read.delim("~/Dropbox/X/chrInf_hg19.txt", header = F, stringsAsFactors = F)
chrlen <- cumsum(as.numeric(chr$V2)) * 1e-6  # mb
names(chrlen) <- c(1:22, "X", "Y")
chrmid <- diff(c(0, chrlen)) * 0.5 + c(0, chrlen[-length(chrlen)])

save(hg19, chrlen, chrmid, file = "Shiny/hg19.rdt")

setwd(hpc)
load("GLM/glmSNP.rdt")
load("GLM/glmIND.rdt")

glm <- do.call(rbind, glm.fit)
glm <- glm[glm$LOD > 6, ]

glm$TYPE <- "SNP"
var.df <- glm[, c("TYPE", "UID", "ID", "CHR", "POS", "LOD")]

glm$TYPE <- "IND"
var.df <- rbind(var.df, glm[, c("TYPE", "UID", "ID", "CHR", "POS", "LOD")])

var.df$POS <- var.df$POS * 1e-3
rownames(var.df) <- NULL

setwd(github)
save(var.df, file = "Shiny/lod.rdt")

# --- TABLE ---
setwd(github)
load("data/tableSNP.rdt")
table0 <- table
load("data/tableIND.rdt")
table <- rbind(table0, table)
save(table, file = "Shiny/table.rdt")

# --- GWAS MANHATTAN ---
glm <- var.df[var.df$LOD > 10, ]
glm$POS <- glm$POS * 1e-3  # mb
glm$CHR1 <- gsub("X", 23, glm$CHR)
glm$CHR1 <- as.numeric(gsub("Y", 24, glm$CHR1))
glm$POS1 <- c(0, chrlen)[glm$CHR1] + glm$POS
glm$COL <- rep("1", nrow(glm))
glm$COL[glm$CHR1 %% 2 == 1] <- "2"

COL <- rep("dodgerblue3", nrow(glm))
COL[glm$CHR1 %% 2 == 1] <- "firebrick1"

setwd(github)
ggplot(glm, aes(x = POS1, y = LOD, color = as.factor(COL))) +
  geom_point(alpha = 0.7) + guides(color = F) +
  theme_bw() + xlab("") + ylab("LOD") +
  scale_color_manual(values = mycol) +
  scale_x_continuous(breaks = chrmid, labels = names(chrlen)) +
  theme(axis.text.x = element_text(color = rep(rev(mycol))))

glm %>%
  ggvis(~POS1, ~LOD, fill = ~COL) %>%
  layer_points() %>%
  add_axis("x", title = "", values = chrmid)

#  tabPanel("Pathway",
#    fluidPage(
#    p(strong("Figure 2:"), "Interaction between transcription factors returned from GWAS and AD-related genes from ALZ forum."),
#    hr(),
#    img(src = "gene_TF.png", height = 512, width = 512)
#    )
#  ),

