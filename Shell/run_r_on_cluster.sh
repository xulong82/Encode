#!/bin/bash
#PBS -l nodes=1:ppn=20,walltime=59:59:59

module load R

# # w/ parameter
# Rscript --no-restore --quiet ~/Dropbox/Network/R1/analysis3.R ${index}

# w/o parameter
# Rscript --no-restore --quiet ~/Dropbox/Network/R1/esb4Df.R
# Rscript --no-restore --quiet ~/Dropbox/Network/R1/liverDf.R
# Rscript --no-restore --quiet ~/Dropbox/Network/R1/melDf.R
  Rscript --no-restore --quiet ~/Dropbox/Network/R1/characterization.R
# Rscript --no-restore --quiet ~/Dropbox/Network/R1/analysis1.R

