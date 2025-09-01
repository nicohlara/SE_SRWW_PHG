#!/bin/sh
#SBATCH --qos=normal
#SBATCH -p atlas
#SBATCH -A guedira_seq_map
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -t 24:00:00
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o "stdout.%x.%j.%N"
#SBATCH -e "stderr.%x.%j.%N"

#genome_path=/project/guedira_seq_map/mwillman/Genome_Assembly_MW/
genome_path=/reference/data/IWGSC/RefSeq_Assemblies/v2.1/
new_path=/90daydata/guedira_seq_map/nico/pangenome/input_data/
#genome='Taes_Hilliard_1.2.fasta'
#genome='Taes_AGS2000_1.1.fasta'
genome=iwgsc_refseqv2.1_assembly.fa
#chromosome='chr1B'
##IWGSC uses different chrom naming
chromosome='Chr1B'
#output='HILLIARD_1B'
#output='AGS2000_1B'
output='CS2.1_1B'

cp $genome_path/$genome $new_path/$genome

module load samtools
samtools faidx ${new_path}/$genome
samtools faidx ${new_path}/$genome ${chromosome} > ${new_path}/${output}.fa
