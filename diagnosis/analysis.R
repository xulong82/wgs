#  Bayesian inference on ADSP project
#  -------------------------------------------------------------------
#  Random effect: Genotype-based Relatedness

library(ggplot2)

cat("--- Extract WAIC data \n")
fname1 <- paste("waic", 1:20, sep = "")
fname2 <- paste(fname1, "rdt", sep = ".")
fname3 <- paste("~/Dropbox/ADSP/Stan", fname2, sep = "/")
load(fname3[1])
myWaic <- get(fname1[1])
for (i in 2:20) {
  if (file.exists(fname3[i])) {
    load(fname3[i])
    myWaic <- c(myWaic, get(fname1[i]))
  }
}

load("/data/xwang/ADSP/R/dbSNP137.rdt")

table(names(myWaic) %in% dbSNP137$ID)  # All TRUE or problem

snp.inf <- dbSNP137[match(names(myWaic), dbSNP137$ID), ]
myWaic <- cbind(snp.inf, myWaic)

save(myWaic, file = "~/Dropbox/ADSP/Stan/myWaic.rdt")

load("/data/xwang/ADSP/R/sscan_clm_rs.rdt")
data <- cbind(sscan[match(myWaic$ID, sscan$ID), ], waic = myWaic[, -c(1:3)])

