#!/bin/bash
# Author: XuLong Wang (xulong.wang@jax.org)

# --- Create compact VCF of the dbSNP137 database from raw

cd /data/xwang/ADSP/dbSNP
file=/data/shared/cga_reference_data/dbsnp_137.hg19.vcf

echo $file

# meta lines
sed '/^##/d' ${file} > temp1

# take the snp id and genotype value columns
awk '{print $3"\t"$1"\t"$2"\t"$4"\t"$5}' temp1 > dbSNP137Hg19.vcf

# remove temporary files
rm temp1

