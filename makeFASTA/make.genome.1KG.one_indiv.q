#!/bin/bash
#$ -m a
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/


-- Note: this slurm script is launched by the script "launch.make.genomes.slurm.sh"

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"
echo "Indiv: ${INDIV}"



ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)


CHR=${ARRY[$SGE_TASK_ID-1]}
echo "chr.q: ${CHR}"

-- path to config.ini file that specifies variables and parameters needed here.
CONFIGPATH='/wynton/group/capra/projects/modern_human_3Dgenome/bin/activeNotebooks/config_makeFASTA_1KG_SNVs.ini'


VCFPTH=$(cat $CONFIGPATH | getSetting 'PATH' 'VCF_PATH')

echo $VCFPTH
-- PFX='ALL.chr'
-- SFX='.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'
--
-- CHRPATH='/wynton/group/capra/data/wynton_databases/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes'
-- VCFHEAD='/wynton/home/capra/egilbertson/data/vcf.header.txt'
-- REFCHRDIR='/wynton/home/capra/egilbertson/data/human_genome/chrms'
--
-- cd /wynton/group/capra/data/wynton_databases/1000_genomes/release/20130502
-- DIRNAME=$(grep ${INDIV} integrated_call_samples_v3.20130502.ALL.panel | awk '{print $3"_"$2"_"$4"_"$1}')
-- echo "Directory: ${DIRNAME}"
--
-- WRKDIR='/wynton/group/capra/projects/modern_human_3Dgenome/data/genomes'
-- cd ${WRKDIR}
--
-- mkdir -p "${DIRNAME}"
-- cd ${DIRNAME}
--
--
-- bash /wynton/group/capra/projects/modern_human_3Dgenome/bin/makeFASTA/make.genome.1KG.one_chromosome.sh ${CHR} ${INDIV} ${VCFPTH} ${PFX} ${SFX} ${DIRNAME} ${CHRPATH} ${WRKDIR} ${VCFHEAD} ${REFCHRDIR}
