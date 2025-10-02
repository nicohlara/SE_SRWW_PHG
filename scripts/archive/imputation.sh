#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=7-00:00:00
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="imputesmall_from_fasta"
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

to_impute=keyfiles/exome_impute.txt

export _JAVA_OPTIONS="-Xmx350G"

module load miniconda3
conda activate /home/nicolas.lara/.conda/envs/phgv2-conda

cd /90daydata/guedira_seq_map/nico/pangenome_multichrom

#phg=phg/bin/phg
phg=../phgv2_v2.4/bin/phg


${phg} list-samples \
	--db-path vcf_dbs \
	--data-set hvcf \
	--output-file output/sample_names_hvcf.txt

hvcf_dir=output/hvcf_export
mkdir ${hvcf_dir}

${phg} export-vcf \
	--db-path vcf_dbs \
	--dataset-type hvcf \
	--sample-file output/sample_names_hvcf.txt \
	--output-dir ${hvcf_dir}

index_dir=output/pangenome_index
index_prefix=soft7_index
#mkdir ${index_dir}

#${phg} rope-bwt-index \
#	--db-path vcf_dbs \
#	--hvcf-dir ${hvcf_dir} \
#	--output-dir ${index_dir} \
#	--index-file-prefix ${index_prefix} \
#	--threads 8

mapped_dir=output/mapped_reads
#mkdir ${mapped_dir}

#${phg} map-reads \
#	--hvcf-dir ${hvcf_dir} \
#	--index ${index_dir}/${index_prefix}.fmd \
#	--key-file ${to_impute} \
#	--output-dir ${mapped_dir}

imputed_hvcf=output/imputed_hvcf
#mkdir ${imputed_hvcf}

#${phg} find-paths \
#	--path-keyfile ${mapped_dir}/pathKeyFile.txt \
#	--hvcf-dir ${hvcf_dir} \
#	--reference-genome output/updated_assemblies/Ref.fa \
#	--path-type haploid \
#	--output-dir ${imputed_hvcf}

imputed_maf=output/imputed_hvcf_maf
#mkdir ${imputed_maf}

#${phg} create-maf-vcf \
#	--db-path vcf_dbs \
#	--bed output/all_ranges.bed \
#	--reference-file output/updated_assemblies/Ref.fa \
#	--maf-dir ${imputed_maf} \
#	-o ${imputed_hvcf}

#${phg} load-vcf \
#	--vcf-dir ${imputed_hvcf} \
#	--db-path vcf_dbs \
#	--threads 10


### work done starting 8-26
h_index=output/index_files
#mkdir ${h_index}

#${phg} rope-bwt-index \
#	--db-path vcf_dbs \
#	--hvcf-dir ${imputed_hvcf} \
#	--output-dir ${h_index} \
#	--index-file-prefix hindex

read_mapping=output/read_mappings
#mkdir ${read_mapping}

#${phg} map-reads \
#    --hvcf-dir output/vcf_files \
#    --index ${h_index}/hindex.fmd \
#    --key-file keyfiles/read_mapping_data.txt \
#    --output-dir ${read_mapping}

imputed_vcf=output/vcf_files_imputed
#mkdir ${imputed_vcf}

#${phg} find-paths \
#    --path-keyfile output/read_mappings/pathKeyFile.txt \
#    --hvcf-dir ${imputed_hvcf} \
#    --reference-genome output/pangenome_index/pangenome.fa \
#    --path-type haploid \
#    --output-dir ${imputed_vcf}

#	--path-keyfile keyfiles/path_finding_data.txt
#	--reference-genome output/updated_assemblies/Ref.fa


## work done pre 8-26
#imputed_snp=output/imputed_snp
#mkdir ${imputed_snp}

#${phg} hvcf2gvcf \
#	--hvcf-dir ${imputed_hvcf} \
#	--db-path vcf_dbs \
#	--reference-file output/updated_assemblies/Ref.fa \
#	--output-dir ${imputed_snp}
