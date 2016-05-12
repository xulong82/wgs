#!/bin/bash
#PBS -l mem=64gb,nodes=1:ppn=20,walltime=72:00:00

module load R

# w/o parameter
# Rscript --no-restore --quiet /data/xwang/bgwas/geno.R
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/glmm/R/sampling_Apoe.R

# w/ parameter
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/glmm/R/optimizing.R ${index}
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/glmm/R/sampling.R ${index}
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/glmm/R/sampling_prior.R ${index}
  Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/glmm/R/sampling_prior_patch.R ${index}
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/glmm/R/optimizing_noK.R ${index}
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/glmm/R/optimizing_prior.R ${index}

