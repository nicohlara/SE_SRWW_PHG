#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=7-00:00:00
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="sam_to_plink"
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


module load samtools
module load bcftools
module load picard
module load plink2
module load conda
conda activate phgv2-conda

dir=/90daydata/guedira_seq_map/nico/pangenome_multichrom/

cd $dir 

input_files=output/alignment_files/*sam

#(printf "@HD\tVN:1.6\tSO:unsorted\n"; samtools dict -H output/updated_assemblies/Ref.fa) > header.sam

#for file in ${input_files}; do
#	name=$(basename $file .sam)
#	cat header.sam $file > fixed_input.sam
#	samtools view -Sb fixed_input.sam | samtools sort -o output_bam/$name.sorted.bam
#done
#rm fixed_input.sam
#rm header.sam

#bam_files=output_bam/*bam
#bcftools mpileup -f output/updated_assemblies/Ref.fa ${bam_files} | bcftools call -mv -Oz -o assemblies.vcf.gz
#bcftools index -f assemblies.vcf.gz

#bcftools norm -f output/updated_assemblies/Ref.fa -Oz -o assemblies_norm.vcf.gz assemblies.vcf.gz
#bcftools index -f assemblies_norm.vcf.gz

#bcftools view -m2 -M2 -v snps assemblies_norm.vcf.gz -Oz -o assemblies_norm_filter.vcf.gz
#bcftools view -s "^output_bam/Ref.sorted.bam" -Oz -o assemblies_norm_filter2.vcf.gz assemblies_norm_filter.vcf.gz
#mv assemblies_norm_filter2.vcf.gz assemblies_norm_filter.vcf.gz
#bcftools index -f assemblies_norm_filter2.vcf.gz


cd haploblocks

#plink2 --vcf assemblies_norm_filter2.vcf.gz --make-bed --allow-extra-chr --out assemblies_plink

##prune the marker set to keep only those markers in approximately perfect linkage equilibrium
##this should in theory give us recombination-based haploblocks to look at
##Really should be run without --bad-ld with a founder population of >50--rerun with the full exome panel at some point?

plink2 --set-missing-var-ids "@_#_cs2" --bfile assemblies_plink --make-bed --allow-extra-chr --out assemblies_plink2
plink2 --bfile assemblies_plink2 --indep-pairwise 1000000 1 0.9 --allow-extra-chr --bad-ld --out recombination_thin

