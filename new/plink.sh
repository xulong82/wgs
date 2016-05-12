#!/bin/bash

cd /data/xwang/adsp3/plink

#------- Note: tped to bed

~/plink_linux_x86_64/plink --tped wgs.tped -tfam wgs.tfam --make-bed --out wgs

#------- Note: sample clustering

~/plink_linux_x86_64/plink --bfile wgs --genome  --out ibs
~/plink_linux_x86_64/plink --bfile wgs --read-genome ibs.genome --out mds --cluster --mds-plot 4

#------- Note: kinship by plink

~/plink_linux_x86_64/plink --bfile wgs --cluster --matrix --out wgs_kin

#------- Note: PCA by plink

~/plink_linux_x86_64/plink --bfile wes --pca --out ${loc}/wes_pca

#------- Note: GWAS

~/plink_linux_x86_64/plink \
  --bfile wes \
  --logistic --out ${loc}/plink \
  --covar ${loc}/wes_cov.txt \
  --covar-number 1-5

# remove "Enriched"

~/plink_linux_x86_64/plink \
  --bfile wes \
  --logistic --out ${loc}/plink_remove_enriched \
  --remove-fam ${loc}/enriched.txt \
  --covar ${loc}/wes_cov.txt \
  --covar-number 1-2

~/plink_linux_x86_64/plink \
  --bfile autosome \
  --logistic --out "$out"/plink_broad \
  --keep-fam /home/xwang/Dropbox/GitHub/wes/docs/remove.txt \
  --covar wes_cov.txt \
  --covar-number 1-2

