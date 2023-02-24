#!/bin/sh
#SBATCH --time=0:01:00
###SBATCH --time=1:00:00

/package/R-base/4.2.0/lib64/R/bin/Rscript  --vanilla reanalyses_DEG_Charly.R

