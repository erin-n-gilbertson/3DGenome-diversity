#!/bin/bash
#SBATCH --job-name=runAkitaComps
#SBATCH --mail-user=erin.gilbertson@vanderbilt.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=0-3:00:00
#SBATCH --nodes=1
#SBATCH --mem=300G #80G
#SBATCH --begin=now
#SBATCH --output=stdout/runAkitaComps_%a.out
#SBATCH --array=1-190%20


echo SBATCH_JOB_NAME: $SBATCH_JOB_NAME
echo SLURM_JOBID:  $SLURM_JOBID

source activate akita

indivs=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'NR==var' listOfPairwiseComps.txt)

python run3dComparisons.DCR.ENG.py "$indivs" > runAkitaComps_"$SLURM_ARRAY_TASK_ID".python.out
