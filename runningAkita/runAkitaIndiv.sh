#!/bin/bash
#SBATCH --job-name=runAkitaIndiv
#SBATCH --mail-user=erin.gilbertson@vanderbilt.edu
#SBATCH --mail-type=FAIL
#SBATCH --time=0-55:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --begin=now
#SBATCH --output=runAkitaOneIndiv_%a.out
#SBATCH --array=1-89%10
#SBATCH --exclude=cn1345,cn1399,cn1426
#!/usr/bin/env python

echo SBATCH_JOB_NAME: $SBATCH_JOB_NAME 
echo SLURM_JOBID:  $SLURM_JOBID

source activate akita

#Identity individual to run Akita on using the listOfIndivs.txt file and array taskid
indiv=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'NR==var' listOfIndivs.txt)

#python runAkitaOnOneIndiv_noHarmonization.py "$indiv"
python runAkitaOnOneIndiv_noHarmonization.py "$indiv"
