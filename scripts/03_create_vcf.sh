#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=1-00:00:00
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="phg_create_vcf"
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

conda init bash
conda activate phgv2-conda
export _JAVA_OPTIONS="-Xmx350G"

cd /90daydata/guedira_seq_map/nico/plinkhaplo_phg

phg=../phgv2_v2.4/bin/phg

#$phg create-ref-vcf \
#	--bed haploblocks/haploblocks.bed \
#        --reference-file output/updated_assemblies/Ref.fa \
#	--reference-name Ref \
#	--db-path vcf_dbs

#mdkir output/vcf_files

#$phg create-maf-vcf \
#    --db-path vcf_dbs \
#    --bed haploblocks/haploblocks.bed \
#    --reference-file output/updated_assemblies/Ref.fa \
#    --maf-dir output/alignment_files \
#    -o output/vcf_files

$phg load-vcf \
	--vcf-dir output/vcf_files \
	--db-path vcf_dbs \
	--threads 10
