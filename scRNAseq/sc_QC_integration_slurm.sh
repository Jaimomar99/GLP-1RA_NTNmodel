#!/bin/bash

#SBATCH --job-name sc_QC_integration
#SBATCH --partition compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=64G
#SBATCH -o /./logs/sc_QC_integration.%A_%a.out

PATH_folder = 

#source conda env
source ""
conda activate ""
DATE=`date +"%Y-%m-%d_%b"`

Rscript ./bin/revised/sc_QC_integration.R\
    > ./logs/sc_QC_integration_${DATE}.log 2>&1