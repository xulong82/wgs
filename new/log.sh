#!/bin/bash

cd ~/Dropbox/GitHub/Adsp/new/log1

files=`ls ./adsp.sh.e*`

for file1 in $files; do
  sed -n '6p' $file1 > ${file1}.x1
  sed -n '15p' $file1 > ${file1}.x2
  paste ${file1}.x1 ${file1}.x2 > ${file1}.x
  rm ${file1}.x[12]
done

files=`ls ./adsp.sh.*.x`

for file1 in $files; do
  cat ${file1} >> log1.txt
done

rm $files

cd ~/Dropbox/GitHub/wes/genotype/log2

files=`ls ./adsp*.sh.e*`

for file1 in $files; do
  sed -n '6p' $file1 > ${file1}.x1
  sed -n '13p' $file1 > ${file1}.x2
  paste ${file1}.x1 ${file1}.x2 > ${file1}.x
  rm ${file1}.x[12]
done

files=`ls ./adsp*.sh.*.x`

for file1 in $files; do
  cat ${file1} >> log2.txt
done

rm $files


