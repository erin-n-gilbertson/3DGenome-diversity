#!/bin/bash
#$ -N runAkitaComps
#$ -l h_rt=72:00:00 #3
#$ -l mem_free=60G #80G300
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/
#$ -t 1-392


echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

source ~/.bash_profile
source ~/.bashrc
#source /wynton/home/capra/egilbertson/envs/akita/bin/activate

load_conda
conda activate modern3d

# change directories
cd /wynton/group/capra/projects/modern_human_3Dgenome/data/pairwise/divergent_windows 

# assign variables
window_list="/wynton/group/capra/projects/modern_human_3Dgenome/data/divergent_windows_exp_distributions.txt"
chr=$(awk -F"," -v row=$SGE_TASK_ID 'NR == row+1 {print $1}' "$window_list")
wndw=$(awk -F"," -v row=$SGE_TASK_ID 'NR == row+1 {print $2}' "$window_list")
script="/wynton/group/capra/projects/modern_human_3Dgenome/bin/comparingAkitaPreds/pairwise_comparisons_by_window_with_dicts.py"


python3 "$script" --chromosome "$chr" --window "$wndw" > /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/runAkitaComps_"$SGE_TASK_ID".python.out
echo "$chr: $wndw comparisons complete"
