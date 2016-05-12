#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=5:00:00

module load parallel
module load vcf-tools
module load bcftools
module load tabix

#--- Note: extract "BOTH_passed" variants: 30 mins

# cd /data/xwang/adsp3/adsp
# dir="/sdata/ADSP/dbGaP-6143/files/vcf_snv/c_all/combined"
 
# file=`find ${dir}/*.chr${chr}.*.vcf.gz`
# vcftools --gzvcf ${file} --keep-filtered "BOTH_passed" --recode --out chr${chr}_both_passed.vcf
  
#--- Note: check sample order 

# cd /data/xwang/adsp3/adsp
# ls chr[1-9][0-9]_*.vcf | parallel "bcftools query -l {} > {}.sp"
# ls *.sp | parallel diff chr10_both_passed.vcf.recode.vcf.sp {}
# rm *.sp

#--- Note: filter variants

# cd /data/xwang/adsp3/golden
# file=chr${chr}_both_passed.vcf.recode.vcf

# vcftools --vcf ../adsp/$file --maf 0.01 --recode --out chr${chr}

#--- Note: concatenate vcf files

# files=`ls *.recode.vcf`
# vcf-concat $files > ../wgs.vcf

#--- Note: make plink files

cd /data/xwang/adsp3
vcftools --vcf wgs.vcf --plink-tped --out ./plink/wgs 
# rm wgs.vcf
 
