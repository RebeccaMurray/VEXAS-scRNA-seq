#!/bin/bash
#SBATCH --job-name=WNN_clustering
#SBATCH --mail-type=all
#SBATCH --mail-user=rmurray@nygenome.org
#SBATCH --partition=pe2
#SBATCH --mem=20g
#SBATCH --output=logs/%x_%j.log

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/rmurray/miniconda3/envs/myR4_RMM_3_17_23

/gpfs/commons/home/rmurray/miniconda3/envs/myR4_RMM_3_17_23/bin/Rscript /gpfs/commons/home/rmurray/rscripts/VEXAS_normal_comparison_MSKCC/R/sample_integration/individual_WNN_clustering.R $1
