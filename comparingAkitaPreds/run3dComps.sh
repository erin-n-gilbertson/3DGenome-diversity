#!/bin/bash
#$ -N runAkitaComps
#$ -l h_rt=1:30:00 #3
#$ -l mem_free=20G #80G300
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/
#$ -t 1-3133756


echo "taskid: ${SGE_TASK_ID}"
echo "JOB_NAME: ${JOB_NAME}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d

source ~/bin/bash_utils/ini_parse

COMP_LIST=$(cat ${CONFIGPATH} | getSetting 'FILE' 'comp_list')



indivs=$(awk -v var="$SGE_TASK_ID" 'NR==var' ${COMP_LIST})

echo $indivs
PYSCRIPT=$(cat ${CONFIGPATH} | getSetting 'BIN' 'run_comps_indiv')


python ${PYSCRIPT} "$indivs" ${CONFIGPATH}> /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/runAkitaComps_"$SGE_TASK_ID".python.out
