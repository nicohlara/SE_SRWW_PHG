

phg initdb --db-path vcf_dbs --gvcf-anchor-gap 10000000 --hvcf-anchor-gap 10000

mkdir output
mkdir output/updated_assemblies

phg prepare-assemblies --keyfile input_data/annotation_keyfile.txt --threads 10 --output-dir output/updated_assemblies

phg agc-compress \
	--db-path vcf_dbs \
	--fasta-list input_data/assemblies_list.txt \
	--reference-file output/updated_assemblies/CS.fa

phg create-ranges --gff input_data/cs2.1_annotation1B_HC.gff3 \
	--reference-file output/updated_assemblies/CS.fa \
	--boundary gene \
	--pad 1000 \
	--range-min-size 500 \
	-o output/ref_ranges.bed

phg align-assemblies \
	--gff input_data/cs2.1_annotation1B_HC.gff3 \
	--reference-file output/updated_assemblies/CS.fa \
	--assembly-file-list input_data/assemblies_list.txt \
	--total-threads 20 \
	--in-parallel 2 \
	--just-ref-prep \
	-o output/alignment_files

phg prepare-slurm-align-file \
    --phg-location /90daydata/guedira_seq_map/nico/pangenome/phg/bin/phg \
    --gff input_data/cs2.1_annotation1B_HC.gff3 \
    --reference-file output/updated_assemblies/CS.fa \
    --reference-sam output/alignment_files/CS.sam \
    --reference-cds-fasta output/alignment_files/ref.cds.fasta \
    --assemblies input_data/assemblies_list.txt \
    --ref-max-align-cov 1 \
    --query-max-align-cov 1 \
    --total-threads 20 \
    --output-dir output/align-assemblies \
    --slurm-command-file output/slurm_align_file.txt \
    -o output/alignment_files

##changed this in prepare-slurm-align-file call
#mv output/alignment_files/ref.cds.fasta output/alignment_files/CD.cds.fasta
##run phg_test_align.sh

phg create-ref-vcf \
	--bed output/ref_ranges.bed \
        --reference-file output/updated_assemblies/CS.fa \
	--reference-name CS \
	--db-path vcf_dbs

phg create-maf-vcf \
    --db-path vcf_dbs \
    --bed output/ref_ranges.bed \
    --reference-file output/updated_assemblies/CS.fa \
    --maf-dir output/alignment_files \
    -o output/vcf_files

phg load-vcf --vcf-dir output/vcf_files \
	--db-path vcf_dbs \
	--threads 10
