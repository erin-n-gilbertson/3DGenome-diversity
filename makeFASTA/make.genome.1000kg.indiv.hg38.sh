# Usage:
# test.sh <CHROMSOME NAME> <INDIVIDUAL IDENTIFIER> <PATH TO DIRECTORY HOLDING VCF FILE>
# Example:
# make.genome.1000.kg.indiv.sh 22 NA19159 '/gpfs51/dors2/capra_lab/data/1000_genomes_project/phase3/vcf/' 'ALL.chr' '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

# /dors/capra_lab/users/rinkerd/projects/3DNeand/bin/make.genome.1000kg.indiv.sh
#cd /dors/capra_lab/users/erin/RotationProject_Akita/data/genomes/1KG/EAS_CHB_male_NA18624
#cd /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/data/genomes
module load Sali CBI gcc gatk bedtools2 samtools htslib

echo "modules loaded"


pwd
#!/bin/bash
### make bed files for all chromosomes lengths

awk '{print $1, "0", $2}' OFS='\t' /wynton/group/capra/data/wynton_databases/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes > hg38.chrom.bed

CHR=$1	# '22'
INDIV=$2	# 'NA19159'
VCFPATH=$3	# '/gpfs51/dors2/capra_lab/data/1000_genomes_project/phase3/hg38_vcf/'
VCFPREFIX=$4	# 'ALL.chr'
VCFSUFFIX=$5 # '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'

INPUTVCF=${VCFPATH}${VCFPREFIX}${CHR}${VCFSUFFIX}

echo "variables defined"
### Find column index for individual

rm tmp${INDIV}.${CHR}
zcat $INPUTVCF | grep -m 1 $INDIV > tmp${INDIV}.${CHR}

if ! [ -s tmp${INDIV}.${CHR} ]; then
	echo "Individual "${INDIV}" not found on chromosome "${CHR}"in 1000 genomes vcf"
	exit
fi


IDX=$(sed -n $'1s/\t/\\\n/gp' tmp${INDIV}.${CHR} | grep -nx $INDIV | cut -d: -f1)
rm tmp${INDIV}.${CHR}
### make new vcf file considering only SNPs, and excluding tri-allelic sites

zcat $INPUTVCF | awk -F '\t|:' -v IDV="$IDX" '/^[^#]/ {if($IDV!="0|0" && $IDV!="2|2" && $IDV!="0|2" && $IDV!="2|0" && $8~/VT=SNP/) print "chr"$1, $2, $3, $4, $5, ".", "." ,$IDV}' OFS='\t' > tmp${CHR}${INDIV}.vcf

cat /wynton/home/capra/egilbertson/data/vcf.header.txt tmp${CHR}${INDIV}.vcf > chr${CHR}_${INDIV}.vcf
rm tmp${CHR}${INDIV}.vcf
bgzip -c chr${CHR}_${INDIV}.vcf > chr${CHR}_${INDIV}.vcf.gz
tabix -p vcf chr${CHR}_${INDIV}.vcf.gz
rm chr${CHR}_${INDIV}.vcf

echo "made new vcf"
### build new genome fasta

gatk FastaAlternateReferenceMaker\
 -R /wynton/home/capra/egilbertson/data/human_genome/chrms/chr${CHR}.fa\
 -V chr${CHR}_${INDIV}.vcf.gz\
 -O chr${CHR}_${INDIV}_hg38_full.fa

echo "build fasta genome"
### Fix GATK output's default fasta headers
sed -i "s/>1/>chr$CHR/g" chr${CHR}_${INDIV}_hg38_full.fa
echo "fix fasta headers"
### remove old dict and index file and remake using the corrected fasta header

rm chr${CHR}_${INDIV}_hg38_full.dict
rm chr${CHR}_${INDIV}_hg38_full.fa.fai

gatk CreateSequenceDictionary -R chr${CHR}_${INDIV}_hg38_full.fa
samtools faidx chr${CHR}_${INDIV}_hg38_full.fa
echo "done"
###
