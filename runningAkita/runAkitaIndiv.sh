#!/bin/bash
#$ -N runAkitaIndiv
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=55:00:00
#$ -l mem_free=10G
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runAkitaIndiv/akita.o
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runAkitaIndiv/akita.e
#$ -t 1
#!/usr/bin/env python

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bashrc
conda activate akita
#Identity individual to run Akita on using the listOfIndivs.txt file and array taskid
indiv=$(awk -v var="$SGE_TASK_ID" 'NR==var' /wynton/group/capra/projects/modern_human_3Dgenome/data/listOfIndivs_all.txt)
echo "Indiv: ${indiv}"
#python runAkitaOnOneIndiv_noHarmonization.py "$indiv"
python /wynton/group/capra/projects/modern_human_3Dgenome/bin/runningAkita/runAkitaOnOneIndiv_noHarmonization.py "$indiv"
