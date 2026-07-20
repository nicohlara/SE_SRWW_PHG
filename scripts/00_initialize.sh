#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=0-12:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="initialize_phg"
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate phgv2-conda

export _JAVA_OPTIONS="-Xmx350G"

# cd /90daydata/guedira_seq_map/nico/plinkhaplo_phg
# phg=../phgv2_v2.4/bin/phg
cd /90daydata/guedira_seq_map/nico/phg_LDblock
#mkdir output


##dir locations
phg=./phg/bin/phg
updated_assemblies=output/updated_assemblies
maf_files=output/maf_files
vcf_files=output/vcf_files

#mkdir ${vcf_files}

##input file locations
assembly_list=data/assemblies_list.txt
gff=/90daydata/guedira_seq_map/nico/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3
bed=/project/guedira_seq_map/nico/pangenome/output/LD_block_algorithm_haploblocks.bed


##initialize a TileDB instance
#${phg} initdb --db-path vcf_dbs
#	--gvcf-anchor-gap 10000000 \
# 	--hvcf-anchor-gap 10000

##update FASTA headers
#${phg} prepare-assemblies \
#	--keyfile data/annotation_keyfile.txt \
#	--threads 10 \
#	--output-dir ${updated_assemblies}

##update chrom names to match throughout
#for fasta in ${updated_assemblies}/*.fa; do
#    sed -i '/^>/ s/chr/Chr/g' "$fasta"
#done

##align assemblies
#${phg} align-assemblies \
#	--gff ${gff} \
#	--reference-file ${updated_assemblies}/Ref.fa \
#	--assembly-file-list ${assembly_list} \
#	-o ${maf_files}

#echo "compress updated assemblies"
#${phg} agc-compress \
# 	--db-path vcf_dbs \
# 	--fasta-list ${assembly_list} \
# 	--reference-file ${updated_assemblies}/Ref.fa

#echo "convert assemblies to vcf"
#${phg} create-ref-vcf \
#	--bed ${bed} \
#	--reference-file ${updated_assemblies}/Ref.fa \
#	--reference-name CS \
#	--db-path vcf_dbs

#echo "get MAF alignments to vcf"
#${phg} create-maf-vcf \
#	--db-path vcf_dbs \
#	--bed ${bed} \
#	--reference-file ${updated_assemblies}/Ref.fa \
#	--maf-dir ${maf_files} \
#	-o ${vcf_files}

###currently breaks on REF
echo "Redoing MAF metrics"
${phg} calc-vcf-metrics \
	--vcf-dir ${vcf_files} \
	--output ${vcf_files}/VCFMetrics_redo.tsv


echo "Load data into DBs"
${phg} load-vcf \
	--vcf-dir ${vcf_files} \
	--db-path vcf_dbs \
	--threads 6

