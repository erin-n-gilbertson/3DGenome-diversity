#!/bin/bash
#$ -N in_silico_mutagenesis_pte_specific_variants
#$ -t 130
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/in_silico_mutagenesis_pte_specific_variants.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/in_silico_mutagenesis_pte_specific_variants.err
#$ -l h_rt=24:00:00
#$ -l mem_free=20G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basenji

# change directories
cd ../../data/in_silico_mutagenesis

# assign variables using the SGE task ID
start=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' ../metadata/splits_for_pte_array.txt)
end=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' ../metadata/splits_for_pte_array.txt)

# run
python3 ../../scripts/3_in_silico_mutagenesis/in_silico_mutagenesis.py --fasta ../reference_fasta/filtered_panTro6.fa --input ../lineage_specific_variants/pte_specific_variants_with_window.txt --start "$start" --end "$end" --out pte_variants_"$SGE_TASK_ID".txt
