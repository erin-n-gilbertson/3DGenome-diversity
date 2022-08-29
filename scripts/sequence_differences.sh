#!/bin/bash
#$ -N sequence_differences
#$ -t 1-3
#$ -M erin.gilbertson@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/sequence_diffs/sequence_differences.out
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/sequence_diffs/sequence_differences.err
#$ -l h_rt=24:00:00
#$ -l mem_free=20G

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
source /wynton/home/capra/egilbertson/envs/akita/bin/activate

# change directories
cd /wynton/group/capra/projects/modern_human_3Dgenome/data/pairwise/sequence

# assign variables using the SGE task ID to get a specific pair of individuals
ind1=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' /wynton/group/capra/projects/modern_human_3Dgenome/data/reference/baseline_sequence_pairs.txt)
ind2=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' /wynton/group/capra/projects/modern_human_3Dgenome/data/reference/baseline_sequence_pairs.txt)

# run
python /wynton/group/capra/projects/modern_human_3Dgenome/bin/scripts/sequence_differences.py --intervals /wynton/group/capra/projects/modern_human_3Dgenome/data/reference/genome_chunks_large.txt --sample_1 /wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/baselines/"$ind1".fa --sample_2 /wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/baselines/"$ind2".fa --sample_1_id "$ind1" --sample_2_id "$ind2"
