#!/bin/sh
#SBATCH --qos=normal
#SBATCH -p atlas
#SBATCH -A guedira_seq_map
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 1:00:00

dir=/90daydata/guedira_seq_map/nico/pangenome/input_data
gff_file=${dir}/iwgsc_refseqv2.1_annotation_200916_HC.gff3
output=${dir}/cs2.1_annotation1B_HC.gff3

grep -E '^(Chr1B)' ${gff_file} > $output
sed -i 's/Chr1B/chr1B/g' $output
