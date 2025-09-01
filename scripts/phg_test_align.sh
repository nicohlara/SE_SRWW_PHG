#!/bin/bash

#SBATCH --account=guedira_seq_map
#SBATCH --time=10:30:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks-per-node=48  # 40 processor core(s) per node X 2 threads per core
#SBATCH --mem=200G  # maximum memory per node
#SBATCH --job-name="10T_anchorwaveV2"
#SBATCH --mail-user=nalara@ncsu.edu  # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load miniconda3

eval "$(conda shell.bash hook)"

conda activate /home/nicolas.lara/.conda/envs/phgv2-conda

echo "All jobs in this array have:"
echo "- SLURM array job id: ${SLURM_ARRAY_JOB_ID}"
echo "- SLURM array task count: ${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM array starting task: ${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM array ending task: ${SLURM_ARRAY_TASK_MAX}"
echo "This job in the array has:"
echo "- SLURM job id: ${SLURM_JOB_ID}"
echo "- SLURM array task id: ${SLURM_ARRAY_TASK_ID}"

INPUTFILE=/90daydata/guedira_seq_map/nico/pangenome_multichrom/output/slurm_align_file.txt

#IFS=$'\n' read -d '' -r -a LINES < ${INPUTFILE}
#LINE=${LINES[$SLURM_ARRAY_TASK_ID]}
while read p; do 
	echo "$p"
	eval ${p}
	if [ $? -eq 0 ]
	  then
	    echo -e "$(date +"%D  %r")\tSuccess: ${LINE}"
	    exit 0
	  else
	    echo -e "$(date +"%D  %r")\tFailed\t${LINE}"
	    echo -e "$(date +"%D  %r")\tJobID\t${SLURM_JOB_ID}"
	    echo -e "$(date +"%D  %r")\tTaskID\t${SLURM_ARRAY_TASK_ID}"
	  exit 1
	fi
done <$INPUTFILE


#eval ${LINE}
#if [ $? -eq 0 ]
#  then
#    echo -e "$(date +"%D  %r")\tSuccess: ${LINE}"
#    exit 0
#  else
#    echo -e "$(date +"%D  %r")\tFailed\t${LINE}"
#    echo -e "$(date +"%D  %r")\tJobID\t${SLURM_JOB_ID}"
#    echo -e "$(date +"%D  %r")\tTaskID\t${SLURM_ARRAY_TASK_ID}"
#  exit 1
#fi
