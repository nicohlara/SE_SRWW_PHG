#!/bin/bash

#SBATCH --account=guedira_seq_map
#SBATCH --time=10:30:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks-per-node=48  # 40 processor core(s) per node X 2 threads per core
#SBATCH --mem=200G  # maximum memory per node
#SBATCH --job-name="10T_anchorwaveV2"
#SBATCH --mail-user=nalara@ncsu.edu  # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load miniconda3

eval "$(conda shell.bash hook)"

conda activate /home/nicolas.lara/.conda/envs/phgv2-conda


cd /90daydata/guedira_seq_map/nico/pangenome_multichrom
#phg/bin/phg align-assemblies \
#        --gff data/iwgsc_refseqv2.1_annotation_200916_HC.gff3 \
#        --reference-file output/updated_assemblies/Ref.fa \
#        --assembly-file-list keyfiles/assemblies_list.txt \
#        --total-threads 18 \
#        --in-parallel 2 \
##        --just-ref-prep \
#        -o output/alignment_files

phg/bin/phg prepare-slurm-align-file \
    --phg-location /90daydata/guedira_seq_map/nico/pangenome/phg/bin/phg \
    --gff data/iwgsc_refseqv2.1_annotation_200916_HC.gff3 \
    --reference-file output/updated_assemblies/Ref.fa \
    --reference-sam output/alignment_files/Ref.sam \
    --reference-cds-fasta output/alignment_files/ref.cds.fasta \
    --assemblies keyfiles/assemblies_list.txt \
    --ref-max-align-cov 1 \
    --query-max-align-cov 1 \
    --total-threads 20 \
    --output-dir output/align-assemblies \
    --slurm-command-file output/slurm_align_file.txt \
    -o output/alignment_files

