#!/bin/bash

cd /data/xwang/ADSP/Plink

plink="$HOME/plink_linux_x86_64/plink"
list="$HOME/Dropbox/GitHub/Adsp/GWAS/list244.txt"
output="$HOME/Dropbox/GitHub/Adsp/GWAS/list244"

$plink --bfile autosome --extract $list --recode-lgen --with-reference --out $output

