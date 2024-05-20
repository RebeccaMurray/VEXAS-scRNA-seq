#!/bin/bash
#SBATCH --job-name=running_rscript
#SBATCH --mail-type=all
#SBATCH --mail-user=rmurray@nygenome.org
#SBATCH --partition=dev
#SBATCH --mem=50g
#SBATCH --cpus-per-task=5
#SBATCH --output=logs/%x_%j.log

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23_seurat_5

module load gdal
module load gsl/2.6
module load zlib/1.2.11 

/gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23_seurat_5/bin/Rscript $1
