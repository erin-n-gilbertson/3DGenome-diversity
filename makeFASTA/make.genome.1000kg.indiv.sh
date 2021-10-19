# Usage:
# test.sh <CHROMSOME NAME> <INDIVIDUAL IDENTIFIER> <PATH TO DIRECTORY HOLDING VCF FILE>
# Example:
# make.genome.1000.kg.indiv.sh 22 NA19159 '/gpfs51/dors2/capra_lab/data/1000_genomes_project/phase3/vcf/' 'ALL.chr' '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

# /dors/capra_lab/users/rinkerd/projects/3DNeand/bin/make.genome.1000kg.indiv.sh
#cd /dors/capra_lab/users/erin/RotationProject_Akita/data/genomes/1KG/EAS_CHB_male_NA18624
cd /dors/capra_lab/users/erin/RotationProject_Akita/data/genomes/1KG/
ml GCC
#ml GATK/4.1.2.0-Java-1.8.0_192
ml GATK
ml BEDTools
ml tabix
ml SAMtools

#!/bin/bash
### make bed files for all chromosomes lengths
cut -f1-2 /gpfs51/dors2/capra_lab/data/dna/human/hg19/hg19.fa.fai > hg19.chrom
awk '{print $1, "0", $2}' OFS='\t' hg19.chrom > hg19.chrom.bed

CHR=$1	# '22'
INDIV=$2	# 'NA19159'
VCFPATH=$3	# '/gpfs51/dors2/capra_lab/data/1000_genomes_project/phase3/vcf/'
VCFPREFIX=$4	# 'ALL.chr'
VCFSUFFIX=$5 # '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

INPUTVCF=${VCFPATH}${VCFPREFIX}${CHR}${VCFSUFFIX}

### Find column index for individual

rm tmp${INDIV}.${CHR}
zcat $INPUTVCF | grep -m 1 $INDIV > tmp${INDIV}.${CHR}

if ! [ -s tmp${INDIV}.${CHR} ]; then
	echo "Individual "${INDIV}" not found on chromosome "${CHR}"in 1000 genomes vcf"
	exit
fi

IDX=$(sed -n $'1s/\t/\\\n/gp' tmp${INDIV}.${CHR} | grep -nx $INDIV | cut -d: -f1)

### make new vcf file considering only SNPs, and excluding tri-allelic sites 

zcat $INPUTVCF | awk -F '\t|:' -v IDV="$IDX" '/^[^#]/ {if($IDV!="0|0" && $IDV!="2|2" && $IDV!="0|2" && $IDV!="2|0" && $8~/VT=SNP/) print "chr"$1, $2, $3, $4, $5, ".", "." ,$IDV}' OFS='\t' > tmp${CHR}${INDIV}.vcf

cat /gpfs51/dors2/capra_lab/users/rinkerd/data/vcf.header.txt tmp${CHR}${INDIV}.vcf > chr${CHR}_${INDIV}.vcf
rm tmp${CHR}${INDIV}.vcf
bgzip -c chr${CHR}_${INDIV}.vcf > chr${CHR}_${INDIV}.vcf.gz
tabix -p vcf chr${CHR}_${INDIV}.vcf.gz
rm chr${CHR}_${INDIV}.vcf

### build new genome fasta

java -jar /gpfs51/dors2/capra_lab/bin/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar FastaAlternateReferenceMaker\
 -R /gpfs51/dors2/capra_lab/data/dna/human/hg19/chr${CHR}.fa\
 -V chr${CHR}_${INDIV}.vcf.gz\
 -O chr${CHR}_${INDIV}_hg19_full.fa

### Fix GATK output's default fasta headers
sed -i "s/>1/>chr$CHR/g" chr${CHR}_${INDIV}_hg19_full.fa

### remove old dict and index file and remake using the corrected fasta header

rm chr${CHR}_${INDIV}_hg19_full.dict
rm chr${CHR}_${INDIV}_hg19_full.fa.fai

java -jar /gpfs51/dors2/capra_lab/bin/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar CreateSequenceDictionary -R chr${CHR}_${INDIV}_hg19_full.fa
samtools faidx chr${CHR}_${INDIV}_hg19_full.fa

### 