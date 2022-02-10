#!/bin/bash
#$ -N runAkitaComps
#$ -l h_rt=1:00:00 #3
#$ -l mem_free=80G #80G300
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/
#$ -t 1


echo "taskid: ${SGE_TASK_ID}"
echo "JOB_NAME: ${JOB_NAME}"

source /wynton/home/capra/egilbertson/envs/akita/bin/activate


indivs=$(awk -v var="$SGE_TASK_ID" 'NR==var' /wynton/group/capra/projects/modern_human_3Dgenome/data/listOfPairwiseComps.txt)

echo $indivs

python /wynton/group/capra/projects/modern_human_3Dgenome/bin/comparingAkitaPreds/run3dComparisons.DCR.ENG.py "$indivs" > runAkitaComps_"$SGE_TASK_ID".python.out
