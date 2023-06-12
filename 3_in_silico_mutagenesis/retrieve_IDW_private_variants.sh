#!/bin/bash
#$ -N retrieve_IDW_private_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/retrieve_IDW_private_variants.out
#$ -e ~/../../../group/capra/projects/pan_3d_genome/scripts/3_in_silico_mutagenesis/retrieve_IDW_private_variants.err
#$ -l h_rt=24:00:00
#$ -l mem_free=30G

# change directories
cd ../../data/IDWs

# set paths
script_path="../../scripts/3_in_silico_mutagenesis/retrieve_private_variants.py"
genotypes_path="../vcfs/one_zero_genotypes.txt"

# run
python3 "$script_path" --genotypes "$genotypes_path" --ids 28 --out chr11_20971520_Jimmie.tmp --region chr11 20971521 22020097
python3 "$script_path" --genotypes "$genotypes_path" --ids 39 --out chr14_26738688_Luky.tmp --region chr14 26738689 27787265
python3 "$script_path" --genotypes "$genotypes_path" --ids 37 --out chr1_72351744_Lara.tmp --region chr1 72351745 73400321
python3 "$script_path" --genotypes "$genotypes_path" --ids 6 --out chr1_168820736_Berta.tmp --region chr1 168820737 169869313
python3 "$script_path" --genotypes "$genotypes_path" --ids 6 --out chr1_169345024_Berta.tmp --region chr1 169345025 170393601
python3 "$script_path" --genotypes "$genotypes_path" --ids 18 --out chr2A_76021760_Coco-chimp.tmp --region chr2A 76021761 77070337
python3 "$script_path" --genotypes "$genotypes_path" --ids 23 --out chr4_82837504_Frederike.tmp --region chr4 82837505 83886081
python3 "$script_path" --genotypes "$genotypes_path" --ids 20 --out chr5_95420416_Desmond.tmp --region chr5 95420417 96468993
python3 "$script_path" --genotypes "$genotypes_path" --ids 9 --out chr6_142606336_Bono.tmp --region chr6 142606337 143654913
python3 "$script_path" --genotypes "$genotypes_path" --ids 3 --out chr7_105906176_Alice.tmp --region chr7 105906177 106954753
python3 "$script_path" --genotypes "$genotypes_path" --ids 5 --out chr8_112197632_Athanga.tmp --region chr8 112197633 113246209
python3 "$script_path" --genotypes "$genotypes_path" --ids 19 --out chr8_128974848_Damian.tmp --region chr8 128974849 130023425

python3 "$script_path" --genotypes "$genotypes_path" --ids 29 53 --out chr10_52428800_Julie-A959_Vaillant.tmp --region chr10 52428801 53477377
python3 "$script_path" --genotypes "$genotypes_path" --ids 27 37 --out chr13_49807360_Ikuru_Lara.tmp --region chr13 49807361 50855937

python3 "$script_path" --genotypes "$genotypes_path" --ids 9 22 33 36 --out chr2A_21495808_Bono_Dzeeta_Kombote_Kumbuka.tmp --region chr2A 21495809 22544385
python3 "$script_path" --genotypes "$genotypes_path" --ids 9 22 33 36 --out chr2A_22020096_Bono_Dzeeta_Kombote_Kumbuka.tmp --region chr2A 22020097 23068673

# concat .tmp files
for f in *.tmp >> IDW_variants_new.txt
