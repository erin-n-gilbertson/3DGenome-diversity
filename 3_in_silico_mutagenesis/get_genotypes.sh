#!/bin/bash
#$ -N get_genotypes
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/get_genotypes.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/get_genotypes.err
#$ -l h_rt=24:00:00
#$ -l mem_free=60G

# load modules
module load CBI
module load bcftools

# change directories
cd ../../data/vcfs

# run
#bcftools merge -0 *_filtered.vcf.gz -O z -o merged_filtered.vcf.gz
#bcftools index -t merged_filtered.vcf.gz
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' merged_filtered.vcf.gz > genotypes.txt
cut --complement -f9,11,12,20,21,25,35,45,50,56,61,62,63,64,68 genotypes.txt > filtered_genotypes.txt
sed 's|0/0|0|g; s|0/1|1|g; s|1/0|1|g; s|1/1|1|g' filtered_genotypes.txt > one_zero_genotypes.txt
