#!/bin/bash
#$ -N map_window_to_pt_specific_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/map_window_to_pt_specific_variants.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/map_window_to_pt_specific_variants.err
#$ -l h_rt=1:00:00
#$ -l mem_free=10G

# load modules
module load CBI
module load bedtools2/2.30.0

# change directories
cd ../../data/lineage_specific_variants

# set variable
pops=('pte' 'pts' 'ptt' 'ptv' )

# generate temp bed files
for p in ${pops[@]}; do awk '{print $1,$2-1,$2,$3,$4}' OFS='\t' "$p"_specific_variants.txt > "$p".bed; done

# intersect
for p in ${pops[@]}; do bedtools intersect -a "$p".bed -b ../metadata/panTro6_windows_with_full_coverage.bed -wa -wb > "$p".tmp; done

# format intersection
for p in ${pops[@]}; do awk '{print $1,$3,$4,$5,$7}' OFS='\t' "$p".tmp > "$p"_specific_variants_with_window.txt; done

# remove other files
for p in ${pops[@]}; do rm "$p".bed && rm "$p".tmp; done
