#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=1-12:00:00
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="impute_from_fasta"
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## to create raw structure for keyfile
#ls input_fasta/*.fast* | awk -F'[/]' '{split($NF,a,"_"); print a[1] "\t" $0}' > keyfile.txt
#make sure to paste in sampleName      filename on top and rename samples if necessary
to_impute=data/GBS_and_exome_impute.txt

export _JAVA_OPTIONS="-Xmx350G"

module load miniconda3
conda activate /home/nicolas.lara/.conda/envs/phgv2-conda

cd /90daydata/guedira_seq_map/nico/plinkhaplo_phg

#phg=phg/bin/phg
phg=../phgv2_v2.4/bin/phg

###imputation prep
#echo "List samples"
#${phg} list-samples \
#	--db-path vcf_dbs \
#	--data-set hvcf \
#	--output-file output/sample_names_hvcf.txt

#echo "Export hvcf files"
hvcf_dir=output/hvcf_export
#mkdir ${hvcf_dir}

#${phg} export-vcf \
#	--db-path vcf_dbs \
#	--dataset-type hvcf \
#	--sample-file output/sample_names_hvcf.txt \
#	--output-dir ${hvcf_dir}

#echo "Index the pangenome for haplotype alignment"
index_dir=output/pangenome_index
index_prefix=soft7_index
#mkdir ${index_dir}

#${phg} rope-bwt-index \
#	--db-path vcf_dbs \
#	--hvcf-dir ${hvcf_dir} \
#	--output-dir ${index_dir} \
#	--index-file-prefix ${index_prefix} \
#	--threads 8

echo "Map reads to pangenome"
read_mapping=output/read_mappings
#mkdir ${read_mapping}

#${phg} map-reads \
#	--hvcf-dir ${hvcf_dir} \
#	--index ${index_dir}/${index_prefix}.fmd \
#	--key-file ${to_impute} \
#	--output-dir ${read_mapping}

imputed_hvcf=output/imputed_hvcf
#mkdir ${imputed_hvcf}

echo "trace haplotype path through graph"
#${phg} find-paths \
#	--path-keyfile ${read_mapping}/pathKeyFile.txt \
#	--hvcf-dir ${hvcf_dir} \
#	--reference-genome output/updated_assemblies/Ref.fa \
#	--path-type haploid \
#	--output-dir ${imputed_hvcf}
#--reference-genome output/pangenome_index/pangenome.fa \ this seems to fail towards end

#echo "index with bcftools"
for vcf in ${imputed_hvcf}/*; do
	bgzip $vcf
	bcftools index ${vcf}.gz
done

echo "Load imputed data into database"
${phg} load-vcf \
	--vcf-dir ${imputed_hvcf} \
	--db-path vcf_dbs \
	--threads 10

echo "Convert haplotype calls to SNP calls"
imputed_snp=output/imputed_snp
#mkdir ${imputed_snp}

${phg} hvcf2gvcf \
	--hvcf-dir ${imputed_hvcf} \
	--db-path vcf_dbs \
	--reference-file output/updated_assemblies/Ref.fa \
	--output-dir ${imputed_snp}
