#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks=48  
#SBATCH --job-name="PHG_load_vcf"
#SBATCH --mail-user=nalara@ncsu.edu  # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
module load miniconda3
module load nextflow
conda activate /home/nicolas.lara/.conda/envs/phgv2-conda

cd /90daydata/guedira_seq_map/nico/pangenome_multichrom

export _JAVA_OPTIONS="-Xmx350G"

##need to change ref_ranges to all_ranges


#./phg/bin/phg agc-compress \
#        --db-path vcf_dbs \
#        --fasta-list keyfiles/assemblies_list.txt \
#        --reference-file output/updated_assemblies/Ref.fa



#./phg/bin/phg align-assemblies \
#        --gff data/iwgsc_refseqv2.1_annotation_200916_HC.gff3 \
#        --reference-file output/updated_assemblies/Ref.fa \
#        --assembly-file-list keyfiles/assemblies_list.txt \
#        --total-threads 18 \
#        --in-parallel 2 \
#        --just-ref-prep \
#        -o output/alignment_files

##
#./phg/bin/phg prepare-slurm-align-file \
#    --phg-location /90daydata/guedira_seq_map/nico/pangenome_multichrom/phg/bin/phg \
#    --gff  data/iwgsc_refseqv2.1_annotation_200916_HC.gff3 \
#    --reference-file output/updated_assemblies/Ref.fa \
#    --reference-sam output/alignment_files/Ref.sam \
#    --reference-cds-fasta output/alignment_files/ref.cds.fasta \
#    --assemblies keyfiles/assemblies_list.txt \
#    --ref-max-align-cov 1 \
#    --query-max-align-cov 1 \
#    --total-threads 20 \
#    --output-dir output/align-assemblies \
#    --slurm-command-file output/slurm_align_file.txt \
#    -o output/alignment_files

#nextflow run /project/guedira_seq_map/nico/pangenome/scripts/align_multiple_assembly.nf

#./phg/bin/phg create-ref-vcf \
#        --bed output/all_ranges.bed \
#        --reference-file output/updated_assemblies/Ref.fa \
#        --reference-name Ref \
#        --db-path vcf_dbs

#mkdir output/vcf_files

#./phg/bin/phg create-maf-vcf \
#    --db-path vcf_dbs \
#    --bed output/all_ranges.bed \
#    --reference-file output/updated_assemblies/Ref.fa \
#    --maf-dir output/alignment_files \
#    -o output/vcf_files

#./phg/bin/phg load-vcf --vcf-dir output/vcf_files \
#        --db-path vcf_dbs \
#        --threads 10
