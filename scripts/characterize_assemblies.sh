#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=1-00:00:00
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="characterize assemblies"
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load seqkit

cd /90daydata/guedira_seq_map/nico2/pangenome_multichrom/output/updated_assemblies/

for assembly in *.fa; do
	echo ${assembly}
	
        seqkit seq -m 400000000 ${assembly} | seqkit stats -
	#seqkit stats ${assembly}
	#awk -v min=50000000 '
	#BEGIN { len=0; name="" }
	#/^>/ {
	#    if (len >= min)
	#        printf "%s\t%d\n", name, len
	#    name = substr($0,2)
	#    len=0
	#    next
	#}
	#{
	#    len += length($0)
	#}
	#END {
    	#if (len >= min)
        #	printf "%s\t%d\n", name, len
	#}' ${assembly}
done > assembly_stats2.txt
