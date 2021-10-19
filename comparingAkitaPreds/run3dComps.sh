#!/bin/bash
#SBATCH --job-name=runAkitaComps
#SBATCH --mail-user=evonne.mcarthur@vanderbilt.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=30G #80G
#SBATCH --begin=now
#SBATCH --output=runAkitaComps_%a.out
#SBATCH --array=1-290%30


echo SBATCH_JOB_NAME: $SBATCH_JOB_NAME 
echo SLURM_JOBID:  $SLURM_JOBID

source activate basenji

indivs=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'NR==var' list_tmp.txt)

python run3dComparisons.DCR.ENG.py "$indivs" > runAkitaComps_"$SLURM_ARRAY_TASK_ID".python.out
