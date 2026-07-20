#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=3-00:00:00  # walltime limit (HH:MM:SS)
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks=48  
#SBATCH --job-name="demultiplex_nextflow"
#SBATCH --mail-user=nalara@ncsu.edu  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


module load miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
#eval "$(command conda shell.bash hook)"
conda activate /home/nicolas.lara/.conda/envs/phgv2-conda

cd /90daydata/guedira_seq_map/nico2

#export _JAVA_OPTIONS="-Xmx350G"

#nextflow scripts/00_vcf_preprocess_combine.nf -resume
#nextflow run /project/guedira_seq_map/nico/pangenome/scripts/00_demultiplex.nf
nextflow run /project/guedira_seq_map/nico/pangenome/scripts/00_align_exome_fastq.nf -resume

squeue -j $SLURM_JOBID
