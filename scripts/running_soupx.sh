#!/bin/bash
#SBATCH --job-name=soupx
#SBATCH --mail-type=all
#SBATCH --mail-user=rmurray@nygenome.org
#SBATCH --partition=pe2
#SBATCH --mem=80g
#SBATCH --output=logs/%x_%j.log
#SBATCH --cpus-per-task=10

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/rmurray/miniconda3/envs/myR4_RMM_3_17_23_soupx

module load gdal
module load gsl/2.6
module load zlib/1.2.11 

/gpfs/commons/home/rmurray/miniconda3/envs/myR4_RMM_3_17_23_soupx/bin/Rscript R/soupx/running_soupx_programmatically.R
