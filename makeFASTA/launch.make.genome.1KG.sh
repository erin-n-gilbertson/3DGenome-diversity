#!/bin/bash
#$ -m a
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/launch.o
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/launch.e
#$ -t 1

# -t 1-2548

# update '-t' flag for however many individuals are in the population list
## the called SLURM script usees an internal array (ARRY) of chromosome names that the slurm array IDs are then used to index
## ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

source ~/bin/bash_utils/ini_parse



###### TODO_Erin
echo $SGE_TASK_ID
path= $(cat ${CONFIGPATH} | getSetting 'PATH' 'LIST_PATH')
list=$(cat ${CONFIGPATH} | getSetting 'PARAMS' 'LIST')
echo "path: $path"
echo "list: $list"


INDIV=$(awk -v var="$SGE_TASK_ID" 'NR==var' $(cat ${CONFIGPATH} | getSetting 'PATH' 'LIST_PATH')/$(cat ${CONFIGPATH} | getSetting 'PARAMS' 'LIST'))
echo $INDIV


qsub -N make.genome.$INDIV -v "INDIV=$INDIV" -l mem_free=80G -t 1-4 -l h_rt=4:00:00 $(cat ${CONFIGPATH} | getSetting 'BIN' 'make_fasta_indiv')
qsub -N make.genome.$INDIV -v "INDIV=$INDIV" -l mem_free=40G -t 5-23 -l h_rt=4:00:00 $(cat ${CONFIGPATH} | getSetting 'BIN' 'make_fasta_indiv')
#

#
