#!/bin/bash
#SBATCH --job-name=running_diffLMM
#SBATCH --mail-type=all
#SBATCH --mail-user=rmurray@nygenome.org
#SBATCH --partition=dev
#SBATCH --mem=150g
#SBATCH --output=logs/%x_%j.log
#SBATCH --cpus-per-task=5

module load anaconda3
source /nfs/sw/anaconda3/anaconda3-10.19/etc/profile.d/conda.sh
conda activate /gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23

/gpfs/commons/home/rmurray/miniconda3/envs/signac_seurat_env_06_04_23/bin/Rscript /gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/R/differential_expression/VEXAS_MUT_vs_WT/lmm/running_lmm_MT_RP_excluded_renormalized_donors_filtered.R
