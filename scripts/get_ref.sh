#!/bin/bash
#$ -N getRef
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -l h_rt=24:00:00
#$ -l mem_free=10G
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/getRef/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/getRef/
#$ -t 1-22


source ~/.bash_profile
source ~/.bashrc
source /wynton/home/capra/egilbertson/envs/akita/bin/activate

module load CBI samtools

cd /wynton/group/capra/data/hg38_ancestral/GAGP
echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"
echo "TASKID: ${SGE_TASK_ID}"

(grep chr$SGE_TASK_ID all_chrs_hg38.txt | head |awk '{printf("%s:%s-%s\n",substr($1,4),$2,$2);}' | while read P; do samtools faidx /wynton/group/capra/data/hg38_fasta/2022-03-14/hg38.fa ${P} ; done) > erin_test/out_chr$SGE_TASK_ID.txt 
