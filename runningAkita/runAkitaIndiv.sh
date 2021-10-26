#!/bin/bash
#$ --N runAkitaIndiv
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt 55:00:00
#$ --l mem_free 10G
#$ -o /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/stdout/runAkitaOneIndiv_%a.out
#$ - 1
#$ -cwd
#!/usr/bin/env python

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"


source activate akita

#Identity individual to run Akita on using the listOfIndivs.txt file and array taskid
indiv=$(awk -v var="$$SGE_TASK_ID" 'NR==var' listOfIndivs_all.txt)
echo "Indiv: ${indiv}"
#python runAkitaOnOneIndiv_noHarmonization.py "$indiv"
python runAkitaOnOneIndiv_noHarmonization.py "$indiv"
