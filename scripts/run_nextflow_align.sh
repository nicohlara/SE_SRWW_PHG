#!/bin/sh
#SBATCH --job-name="nextflow_alignment"
#SBATCH --qos=normal
#SBATCH -p atlas
#SBATCH -A guedira_seq_map
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -t 7-00:00:00
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o "stdout.%x.%j.%N"                                                                                                                                                                                     #SBATCH -e "stderr.%x.%j.%N"

module load nextflow
module load miniconda3

cd /projects/guedira_seq_map/nico/pangenome
nextflow scripts/align_multiple_assembly.nf
