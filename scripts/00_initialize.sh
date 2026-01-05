#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=1-00:00:00
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="initialize_phg"
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

conda activate phgv2-conda

cd /90daydata/guedira_seq_map/nico/plinkhaplo_phg
phg=../phgv2_v2.4/bin/phg

##initialize a TileDB instance
${phg} initdb --db-path vcf_dbs \
	--gvcf-anchor-gap 10000000 \
	--hvcf-anchor-gap 10000

mkdir output
mkdir output/updated_assemblies

##update FASTA headers
#${phg} prepare-assemblies \
#	--keyfile keyfiles data/annotation_keyfile.txt \
#	--threads 10 \
#	--output-dir output/updated_assemblies

##compress updated assemblies
#phg agc-compress \
#	--db-path vcf_dbs \
#	--fasta-list data/assemblies_list.txt \
#	--reference-file output/updated_assemblies/CS.fa

#mkdir output/align_assemblies
##prepare slurm file for parallel alignments via nextflow
#phg prepare-slurm-align-file \
#    --phg-location /90daydata/guedira_seq_map/nico/phgv2_v2.4/bin/phg \
#    --gff input_data/cs2.1_annotation1B_HC.gff3 \
#    --reference-file output/updated_assemblies/CS.fa \
#    --reference-sam output/alignment_files/CS.sam \
#    --reference-cds-fasta output/alignment_files/ref.cds.fasta \
#    --assemblies data/assemblies_list.txt \
#    --ref-max-align-cov 1 \
#    --query-max-align-cov 1 \
#    --total-threads 20 \
#    --output-dir output/align-assemblies \
#    --slurm-command-file output/slurm_align_file.txt \
#    -o output/alignment_files
