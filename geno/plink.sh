#!/bin/bash

# MAKE PLINK TPED FILE
# Author: XuLong Wang (xulong.wang@jax.org)

par=$1

varInf="biInf_chr$par.txt"
geno="biGeno_chr$par.txt"
tped="chr$par.tped"
bed="chr$par"

cd /data/xwang/ADSP/Plink

# nf=`awk 'NR==1 {print NF; exit}' ../biVar/$geno`
# nf=576

awk '{print $1, $5, $6}' ../biVar/$varInf > $par.0
awk '{$8 = 0} {print $2, $1, $8, $3}' ../biVar/$varInf > $tped

for idx in `seq 1 576`; do
  awk -v col=$idx '{print $col}' ../biVar/$geno > $par.1
  paste $par.0 $par.1 > $par.2
  awk '$4==0 {$5=$2; $6=$2} 
       $4==1 {$5=$2; $6=$3} 
       $4==2 {$5=$3; $6=$3} 
       {print $5, $6}' $par.2 > $par.3
  paste $tped $par.3 > $par.4
  mv $par.4 $tped
  rm $par.[123]
done
rm $par.0

