#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks=48  
#SBATCH --job-name="demultiplex_nextflow"
#SBATCH --mail-user=nalara@ncsu.edu  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniconda3
#module load nextflow
conda activate /home/nicolas.lara/.conda/envs/imp_2

cd /90daydata/guedira_seq_map/nico2

#export _JAVA_OPTIONS="-Xmx350G"

#nextflow scripts/00_vcf_preprocess_combine.nf -resume
nextflow run 00_demultiplex.nf

squeue -j $SLURM_JOBID
