#!/bin/bash
#PBS -l nodes=1:ppn=20,walltime=36:00:00

module load R

# w/o parameter
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/Adsp/new/vcf2r.R

# w/ parameter
  Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/Adsp/new/sampling.R ${chr}
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/Adsp/new/optimizing.R ${chr}

