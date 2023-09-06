#!/bin/bash
#$ -N in_silico_mutagenesis__idw_rare_alternate
#$ -M erin.gilbertson@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/in_silico_mutagenesis/in_silico_mutagenesis_idw_rare_alternate.out
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/in_silico_mutagenesis/in_silico_mutagenesis_idw_rare_alternate.err
#$ -l h_rt=24:00:00
#$ -l mem_free=20G

# load conda environment
echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d

# change directories
cd ../../data/in_silico_mutagenesis/IDW


echo "SGE_TASK_ID:  ${SGE_TASK_ID}"
# assign variables using the SGE task ID
start=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' ISM_splits_idw_rare_alternate.txt)
end=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' ISM_splits_idw_rare_alternate.txt)
echo "start: ${start}"
echo "end: ${end}"

# run
python3 ../../../bin/3_in_silico_mutagenesis/in_silico_mutagenesis.py --fasta hg38_reference  --input ../../IDWs/IDW_variants_rare_alternate.txt  --start "$start" --end "$end" --out ism_scores_idw_rare_alternate_"$SGE_TASK_ID".txt
