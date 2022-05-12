#!/bin/bash
#$ -N runAkitaIndiv
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=55:00:00
#$ -l mem_free=10G
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runAkitaIndiv/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runAkitaIndiv/
#$ -t 1

#!/usr/bin/env python
#qsub runAkitaIndiv.sh

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
source /wynton/home/capra/egilbertson/envs/akita/bin/activate

source ~/bin/bash_utils/ini_parse
CONFIGPATH='/wynton/group/capra/projects/modern_human_3Dgenome/bin/activeNotebooks/config_makeFASTA_1KG_SNVs.ini'

INDIV_List=$(cat ${CONFIGPATH} | getSetting 'PATH' 'ind_list')
#Identity individual to run Akita on using the listOfIndivs.txt file and array taskid
indiv=$(awk -v var="$SGE_TASK_ID" 'NR==var' ${INDIV_LIST})
echo "Indiv: ${indiv}"


PYSCRIPT=$(cat ${CONFIGPATH} | getSetting 'BIN' 'run_akita_indiv')

#python PYSCRIPT "$indiv" ${CONFIGPATH}
