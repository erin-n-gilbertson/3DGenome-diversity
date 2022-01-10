#!/bin/bash
#$ -M erin.gilbertson@ucsf.edu
#$ -m a
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/${INDIV}.o
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/${INDIV}.e


## Note: this slurm script is launched by the script "launch.make.genomes.slurm.sh"

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"
echo "Indiv: ${INDIV}"

ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

CHR=${ARRY[$SGE_TASK_ID]}

VCFPTH='/wynton/group/capra/data/wynton_databases/1000_genomes/release/20190312_biallelic_SNV_and_INDEL/'
PFX='ALL.chr'
SFX='.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'

cd /wynton/group/capra/data/wynton_databases/1000_genomes/release/20130502
DIRNAME=$(grep ${INDIV} integrated_call_samples_v3.20130502.ALL.panel | awk '{print $3"_"$2"_"$4"_"$1}')
echo $DIRNAME
cd /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/data/genomes
mkdir -p "${DIRNAME}"
cd ${DIRNAME}


bash /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/bin/makeFASTA/make.genome.1000kg.indiv.hg38.sh ${CHR} ${INDIV} ${VCFPTH} ${PFX} ${SFX} ${DIRNAME}
