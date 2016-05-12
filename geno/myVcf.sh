#!/bin/bash

# Make compact VCF files

cd /data/xwang/ADSP/VCF

files=`find /sdata/ADSP/processed/final/14*/*_filtered.vcf.gz`

for name1 in $files; do
    name2=`basename $name1`
    name3=${name2/_variants_filtered.vcf.gz/}

    echo $name3
    cp $name1 ./

    # unzip the gz
    gzip -d $name2

    # delete meta lines
    sed '/^#/d' ${name2/.gz/} > temp1

    # delete the LowCoverage entries
    grep -v "LowCoverage" temp1 > temp2

    # delete the VeryLowQual entries
    grep -v "VeryLowQual" temp2 > temp3

    # take the snp id and genotype value columns
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10}' temp3 > temp4

    # keep genotype info, discard others
    awk 'gsub(":.*", "", $6) {print $0}' temp4 > temp5

    # header
    awk 'BEGIN {print "CHR POS ID REF ALT GT"} {print $0}' temp5 > ${name3}.vcf

    # remove temporary files
    rm ${name2/.gz}
    rm temp[12345]
done
