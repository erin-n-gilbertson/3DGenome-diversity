#!/bin/bash
#$ -N in_silico_mutagenesis_idw_rare_reference
#$ -M erin.gilbertson@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/in_silico_mutagenesis/in_silico_mutagenesis_idw_all_variants.out
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/in_silico_mutagenesis/in_silico_mutagenesis_idw_all_variants.err
#$ -l h_rt=24:00:00
#$ -l mem_free=20G


# usage qsub in_silico_mutagenesis_idw_all_variants_few_windows.sh -v INPUT=[input_file],OUT="[output_name]"
# load conda environment
echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d

# change directories
cd ../../data/in_silico_mutagenesis/all_variants_few_windows



echo "SGE_TASK_ID:  ${SGE_TASK_ID}"
# assign variables using the SGE task ID
INPUT=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' ISM_variables.txt)
OUT=$(cut -d '_' -f 1,2 <<< "$INPUT")

echo "INPUT: ${INPUT}"
echo "OUT: ${OUT}"

# run
python3 ../../../bin/3_in_silico_mutagenesis/in_silico_mutagenesis.py --fasta hg38_reference  --input ../../IDWs/${INPUT} --out ism_scores_all_variants_${OUT}.txt
