#!/bin/bash
#$ -m a
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/


## Note: this slurm script is launched by the script "launch.make.genomes.slurm.sh"
source ~/bin/bash_utils/ini_parse

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"
echo "Indiv: ${INDIV}"





## path to config.ini file that specifies variables and parameters needed here.
## THIS IS THE ONLY THING THAT SHOULD HAVE TO CHANGE!!!

CONFIGPATH='/wynton/group/capra/projects/modern_human_3Dgenome/bin/activeNotebooks/config_makeFASTA_1KG_SNVs.ini'
echo "configpath: ${CONFIGPATH}"

VCFPTH=$(cat ${CONFIGPATH} | getSetting 'PATH' 'vcf_path')


PFX=$(cat ${CONFIGPATH} | getSetting 'PARAMS' 'vcf_pfx')
SFX=$(cat ${CONFIGPATH} | getSetting 'PARAMS' 'vcf_sfx')

CHRPATH=$(cat ${CONFIGPATH} | getSetting 'PATH' 'chr_path')
VCFHEAD=$(cat ${CONFIGPATH} | getSetting 'PATH' 'vcf_head')
REFCHRDIR=$(cat ${CONFIGPATH} | getSetting 'PATH' 'ref_path')


ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)


CHR=${ARRY[$SGE_TASK_ID-1]}
echo "chr.q: ${CHR}"

cd /wynton/group/capra/data/wynton_databases/1000_genomes/release/20130502
DIRNAME=$(grep ${INDIV} integrated_call_samples_v3.20130502.ALL.panel | awk '{print $3"_"$2"_"$4"_"$1}')
echo "Directory: ${DIRNAME}"

WRKDIR=$(cat ${CONFIGPATH} | getSetting 'PATH' 'fasta_out_path')
cd ${WRKDIR}

mkdir -p "${DIRNAME}"
cd ${DIRNAME}


bash $(cat ${CONFIGPATH} | getSetting 'BIN' 'make_fasta_chrm') ${CHR} ${INDIV} ${VCFPTH} ${PFX} ${SFX} ${DIRNAME} ${CHRPATH} ${WRKDIR} ${VCFHEAD} ${REFCHRDIR}
