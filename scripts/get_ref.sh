#!/bin/bash

source ~/.bash_profile
source ~/.bashrc
source /wynton/home/capra/egilbertson/envs/akita/bin/activate

module load CBI samtools

cd /wynton/group/capra/data/hg38_ancestral/GAGP
(grep -v "#" all_chrs_hg38_strip_chr.txt | awk '{printf("%s:%s-%s\n",$1,$2,$2);}' | while read P; do samtools faidx /wynton/group/capra/data/hg38_fasta/2022-03-14/hg38.fa ${P} ; done) > out.txt 
