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

module load miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate phgv2-conda

export _JAVA_OPTIONS="-Xmx350G"

cd /90daydata/guedira_seq_map/nico/phg_LDblock

##dir locations
phg=./phg/bin/phg
updated_assemblies=output/updated_assemblies
hvcf_dir=output/hvcf_export
index_dir=output/pangenome_index
index_prefix=soft7_index
#mkdir ${hvcf_dir}
#mkdir ${index_dir}

##imputation parameters
to_impute=data/founder_impute_test.txt
project=founder_impute_exome_gbs
read_mapping=output/${project}/read_mappings
imputed_hvcf=output/${project}/imputed_hvcf
imputed_snp=output/${project}/imputed_snp
#mkdir output/${project}
#mkdir ${read_mapping}
#mkdir ${imputed_hvcf}
#mkdir ${imputed_snp}


###imputation prep

##something at the previous step seems a little odd, have to manually move the Ref hvcf out of the directory to proceed.
#echo "Listing samples"
#echo ""
#${phg} list-samples \
#	--db-path vcf_dbs \
#	--data-set hvcf \
#	--output-file output/sample_names_hvcf.txt

#echo "Exporting vcf"
#echo ""
#${phg} export-vcf \
#	--db-path vcf_dbs \
#	--dataset-type hvcf \
#	--sample-file output/sample_names_hvcf.txt \
#	--output-dir ${hvcf_dir}

#echo "Indexing pangenome"
#echo ""
#${phg} rope-bwt-index \
#	--db-path vcf_dbs \
#	--hvcf-dir ${hvcf_dir} \
#	--output-dir ${index_dir} \
#	--index-file-prefix ${index_prefix} \
#	--threads 8


###imputation steps with new data starts here
#echo "Mapping reads"
#echo ""

#${phg} map-reads \
#	--hvcf-dir ${hvcf_dir} \
#	--index ${index_dir}/${index_prefix}.fmd \
#	--key-file ${to_impute} \
#	--min-mem-length 88 \
#	--threads 5 \
#	--output-dir ${read_mapping}

echo "Finding paths"
echo ""

#${phg} find-paths \
#	--path-keyfile ${read_mapping}/pathKeyFile.txt \
#	--hvcf-dir ${hvcf_dir} \
#	--reference-genome output/updated_assemblies/Ref.fa \
#	--path-type haploid \
#	--output-dir ${imputed_hvcf}

for vcf in ${imputed_hvcf}/*.h.vcf; do
	bgzip $vcf
	bcftools index ${vcf}.gz
done

${phg} load-vcf \
	--vcf-dir ${imputed_hvcf} \
	--db-path vcf_dbs \
	--threads 5

#${phg} hvcf2gvcf \
#	--hvcf-dir ${imputed_hvcf} \
#	--db-path vcf_dbs \
#	--reference-file output/updated_assemblies/Ref.fa \
#	--output-dir ${imputed_snp}

