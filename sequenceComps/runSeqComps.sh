#!/bin/bash
#$ -N runSeqComps
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=24:00:00
#$ -l mem_free=10G
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runSeqComps/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runSeqComps/
#$ -t 1-2504


echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
source /wynton/home/capra/egilbertson/envs/akita/bin/activate

source ~/bin/bash_utils/ini_parse

echo "config: ${CONFIGPATH}"

COMP_LIST=$(cat ${CONFIGPATH} | getSetting 'FILE' 'comp_list')

indivs=$(awk -v var="$SGE_TASK_ID" 'NR==var' ${COMP_LIST})
echo $indivs

PYSCRIPT=$(cat ${CONFIGPATH} | getSetting 'BIN' 'run_comps_indiv')

python ${PYSCRIPT} "$indivs" ${CONFIGPATH}> /wynton/group/capra/projects/modern_human_3Dgenome/stdout/runSeqComps/runSeqComps_"$SGE_TASK_ID".python.out
