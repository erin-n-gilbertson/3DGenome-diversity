#!/bin/bash
#$ -M erin.gilbertson@ucsf.edu
#$ -m abe
#$ -o /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/stdout
#$ -e /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/stdout


## Note: this slurm script is launched by the script "launch.make.genomes.slurm.sh"

echo "JOB_NAME: ${JOB_NAME}"
echo "JOBID:  ${JOB_ID}"

ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

CHR=${ARRY[$SGE_TASK_ID]}

cd /wynton/group/capra/data/wynton_databases/1000_genomes/release/20130502
DIRNAME=$(grep ${INDV} integrated_call_samples_v3.20130502.ALL.panel | awk '{print $3"_"$2"_"$4"_"$1}')
mkdir -p ${DIRNAME}
cd ${DIRNAME}


bash /dors/capra_lab/users/erin/RotationProject_Akita/scripts/makeFASTA/make.genome.1000kg.indiv.sh ${CHR} ${INDV} ${VCFPTH} ${PFX} ${SFX}
