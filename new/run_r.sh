#!/bin/bash
#PBS -l nodes=1:ppn=2,walltime=6:00:00

module load R
# module load R/3.2.3

# w/o parameter
  Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/wgs/new/meta2.R
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/Adsp/new/vcf2r.R

# w/ parameter
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/wgs/new/sampling.R ${chr}
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/wgs/new/sampling_prior.R ${chr}
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/wgs/new/optimizing.R ${chr}

