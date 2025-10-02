#!/bin/bash
#SBATCH --account=guedira_seq_map
#SBATCH --time=7-00:00:00
#SBATCH --partition=atlas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="fasta_stats"
#SBATCH --mail-user=nalara@ncsu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

dir=/90daydata/guedira_seq_map/nico/pangenome_multichrom/
input_files=${dir}*fast*/*

# genome size constant
size=15070000000

# print header
echo -e "File\tCoverage"

for file in $input_files; do
    # use zcat if compressed, cat otherwise
    if [[ "$file" == *.gz ]]; then
        reader="zcat"
    else
        reader="cat"
    fi

    # number of reads (FASTQ has 4 lines per read)
    count=$($reader "$file" | wc -l)
    count=$((count / 4))

    # average read length: sequence lines are every 4th line starting at line 2
    avg_len=$($reader "$file" | awk 'NR%4==2 {sum+=length($0); n++} END {if(n>0) print sum/n; else print 0}')

    # coverage
    coverage=$(awk -v c="$count" -v l="$avg_len" -v s="$size" \
        'BEGIN {printf "%.4f", (c*l)/s}')

    # print result
    echo -e "$(basename "$file")\t$count\t${avg_len}\t$size\t$coverage"
done

