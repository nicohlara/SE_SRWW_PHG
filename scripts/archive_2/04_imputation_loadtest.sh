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

to_impute=keyfiles/SunRILs2_impute.txt
project=SunRILs_population


export _JAVA_OPTIONS="-Xmx350G"

module load miniconda3
eval "$(command conda shell.bash hook)"
conda activate /home/nicolas.lara/.conda/envs/phgv2-conda

cd /90daydata/guedira_seq_map/nico2/pangenome_multichrom

#phg=phg/bin/phg
phg=../phgv2_v2.4/bin/phg

imputed_hvcf=output/${project}_imputed_hvcf

imputed_hvcf_loading=${imputed_hvcf}_loading
#mkdir ${imputed_hvcf_loading}
##for testing purposes, pull N lines from imputed_hvcf into temp folder
##updated to only pull lines with BLUEs
N=1000
BLUE=/project/guedira_seq_map/nico/pangenome/data/blues.csv
#find ${imputed_hvcf} -maxdepth 1 -type f -name "*.h.vcf.gz" \
awk -F',' 'NR==1 {
    for (i=1;i<=NF;i++) if ($i=="Entry") col=i
    next
  }
  { print $col }' "${BLUE}" \
  | sed 's/$/.h.vcf.gz/' \
  | while read -r fname; do
      fullpath="${imputed_hvcf}/${fname}"
      [[ -f "$fullpath" ]] && echo "$fullpath"
  done \
  | shuf \
  | head -n "$N" \
  | while read -r vcf; do
      csi="${vcf}.csi"
      mv "$vcf" "$csi" "${imputed_hvcf_loading}/"
    done

##calculate time for N samples with M threads
M=5
start=$SECONDS
${phg} load-vcf \
        --vcf-dir ${imputed_hvcf_loading} \
        --db-path vcf_dbs \
	--threads ${M}
duration=$(( SECONDS - start))
echo "${N} samples loaded in $duration seconds with thread = $M"


imputed_hvcf_loaded=${imputed_hvcf}_loaded
mv ${imputed_hvcf_loading}/* ${imputed_hvcf_loaded}

#imputed_snp=output/${project}_imputed_snp
#mkdir ${imputed_snp}

#${phg} hvcf2gvcf \
#        --hvcf-dir ${imputed_hvcf} \
#        --db-path vcf_dbs \
#        --reference-file output/updated_assemblies/Ref.fa \
#        --output-dir ${imputed_snp}

