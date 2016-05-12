#!/bin/bash

cd /data/xwang/ADSP/Plink

# MAKE AUTOSOME TPED
# rm autosome.tped
# for idx in `seq 1 22`; do
#   chr=chr$idx
#   echo $chr
#   tped="$chr.tped"
#   cat $tped >> autosome.tped
# done
 
# TPED TO BED 
# ~/plink_linux_x86_64/plink --tped autosome.tped -tfam adsp.tfam --make-bed --out autosome

# KING'S KINSHIP
# ~/king -b autosome.bed --homo --prefix ../Kinship/autosome

# LEAVE ONE OUT TPED & KINSHIP
# par=$1
# chr=chr$par
# nochr=no_chr$par
# echo $nochr
# awk -v arg=$par '$1 == arg {next} {print $0}' autosome.tped > $nochr.tped
# ~/plink_linux_x86_64/plink --tped $nochr.tped -tfam adsp.tfam --make-bed --out $nochr
# ~/king -b $nochr.bed --homo --prefix ../Kinship/$nochr

# LEAVE TWO OUT TPED & KINSHIP
  par=$1
  chrA=chr$par
  nochrA=no_chr$par
  for idx in `seq $par 22`; do
    nochrAB=${nochrA}chr$idx
    echo $nochrAB
    awk -v arg=$idx '$1 == arg {next} {print $0}' $nochrA.tped > $nochrAB.tped
    ~/plink_linux_x86_64/plink --tped $nochrAB.tped -tfam adsp.tfam --make-bed --out $nochrAB
    ~/king -b $nochrAB.bed --homo --prefix ../Kinship/$nochrAB
  done

# LD-based SNP pruning
# /home/xwang/plink_linux_x86_64/plink --tped autosome.tped -tfam adsp.tfam --indep 50 5 2

# LD calculation
# ~/plink_linux_x86_64/plink --bfile autosome --ld 1-10250 1-10257
# ~/plink_linux_x86_64/plink --bfile autosome -r2 --ld-snp-list ~/Dropbox/GitHub/Adsp/GWAS/hitList.txt \
#                            --ld-window-kb 1000 \
# 			   --ld-window 10000 \
# 			   --ld-window-r2 0

# awk 'NR==FNR {h[$3] = $7; next} {print $0, h[$2]}' SRR1103888.vcf autosome.map > 1
# awk '$5 {$6 = $5} !$5 {$6 = "NN"} {print $6}'  1 > 2

