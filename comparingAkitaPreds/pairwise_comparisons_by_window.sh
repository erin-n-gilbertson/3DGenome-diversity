#!/bin/bash
#$ -N runAkitaComps
#$ -l h_rt=1:30:00 #3
#$ -l mem_free=60G #80G300
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/comps/
#$ -t 1-392

# change directories
cd /wynton/group/capra/projects/modern_human_3Dgenome/data/pairwise/divergent_windows

# assign variables
window_list="/wynton/group/capra/projects/modern_human_3Dgenome/data/divergent_windows_exp_distributions.txt"
chr=$(awk -F"," -v row=$SGE_TASK_ID 'NR > row {print $1}' "$window_list")
wndw=$(awk -v row=$SGE_TASK_ID 'NR > row {print $2}' "$window_list")
script="/wynton/group/capra/projects/modern_human_3Dgenome/bin/comparingAkitaPreds/pairwise_comparisons_by_windows.py"


python3 "$script" --chromosome "$chr" --window "$window"
echo "$chr: $wndw comparisons complete"
