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

#to_impute=keyfiles/exome_impute.txt
#to_impute=keyfiles/exome_round2_GBS_impute.txt
#to_impute=keyfiles/round2_GBS_impute.txt
to_impute=keyfiles/SunRILs2_impute.txt
project=SunRILs_population


export _JAVA_OPTIONS="-Xmx350G"

module load miniconda3
eval "$(command conda shell.bash hook)"
conda activate /home/nicolas.lara/.conda/envs/phgv2-conda

cd /90daydata/guedira_seq_map/nico2/pangenome_multichrom

#phg=phg/bin/phg
phg=../phgv2_v2.4/bin/phg

###imputation prep
#${phg} list-samples \
#	--db-path vcf_dbs \
#	--data-set hvcf \
#	--output-file output/sample_names_hvcf.txt

hvcf_dir=output/hvcf_export
#mkdir ${hvcf_dir}

#${phg} export-vcf \
#	--db-path vcf_dbs \
#	--dataset-type hvcf \
#	--sample-file output/sample_names_hvcf.txt \
#	--output-dir ${hvcf_dir}

index_dir=output/pangenome_index
index_prefix=soft7_index
#mkdir ${index_dir}

#${phg} rope-bwt-index \
#	--db-path vcf_dbs \
#	--hvcf-dir ${hvcf_dir} \
#	--output-dir ${index_dir} \
#	--index-file-prefix ${index_prefix} \
#	--threads 8


###imputation steps with new data starts here

#read_mapping=output/${project}_read_mappings
#mkdir ${read_mapping}
read_mapping=output/read_mappings_SunRILs_GBS

#${phg} map-reads \
#	--hvcf-dir ${hvcf_dir} \
#	--index ${index_dir}/${index_prefix}.fmd \
#	--key-file ${to_impute} \
#	--min-mem-length 88 \
#	--threads 48 \
#	--output-dir ${read_mapping}

imputed_hvcf=output/${project}_imputed_hvcf
#mkdir ${imputed_hvcf}
#imputed_hvcf=utput/imputed_hvcf_SunRILs_GBS


#${phg} find-paths \
#	--path-keyfile keyfiles/SunRILs2_readmappingKeyfile.tsv \
#	--hvcf-dir ${hvcf_dir} \
#	--reference-genome output/updated_assemblies/Ref.fa \
#	--path-type haploid \
#	--output-dir ${imputed_hvcf}
#        --path-keyfile ${read_mapping}/pathKeyFile.txt \

for vcf in ${imputed_hvcf}/*.h.vcf; do
	bgzip $vcf
	bcftools index ${vcf}.gz
done

${phg} load-vcf \
	--vcf-dir ${imputed_hvcf} \
	--db-path vcf_dbs \
	--threads 48

imputed_snp=output/imputed_snp_SunRILs_GBS
imputed_snp=output/${project}_imputed_snp
mkdir ${imputed_snp}

${phg} hvcf2gvcf \
	--hvcf-dir ${imputed_hvcf} \
	--db-path vcf_dbs \
	--reference-file output/updated_assemblies/Ref.fa \
	--output-dir ${imputed_snp}

