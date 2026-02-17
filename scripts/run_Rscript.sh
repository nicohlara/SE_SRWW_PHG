#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks=48  
#SBATCH --job-name="Run accuracy assessment"
#SBATCH --mail-user=nalara@ncsu.edu  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

##activate env
module load miniconda3
eval "$(command conda shell.bash hook)"
conda activate /home/nicolas.lara/.conda/envs/r-phg


module load openblas
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib/jvm/lib/server:$LD_LIBRARY_PATH

#Rscript /project/guedira_seq_map/nico/pangenome/scripts/05_haplotype_GRM.R
Rscript /project/guedira_seq_map/nico/pangenome/scripts/05_accuracy_assessment.R
#Rscript /project/guedira_seq_map/nico/pangenome/scripts/06_GP_testing.R
#Rscript /project/guedira_seq_map/nico/pangenome/scripts/06_haplotype_GWAS.R
#Rscript /project/guedira_seq_map/nico/pangenome/scripts/07_summary_stats.R
