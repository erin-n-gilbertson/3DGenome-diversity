#!/bin/bash
#SBATCH --mail-user=erin.gilbertson@vanderbilt.edu
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/dors/capra_lab/users/erin/RotationProject_Akita/data/genomes/makeFASTAout/array_job_slurm_%A_%a.out
#SBATCH --reservation=AtRiskUnstableEnvironment

## Note: this slurm script is launched by the script "launch.make.genomes.slurm.sh"

echo SBATCH_JOB_NAME: $SBATCH_JOB_NAME 
echo SLURM_JOBID:  $SLURM_JOBID

ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

CHR=${ARRY[$SLURM_ARRAY_TASK_ID]}

cd /dors/capra_lab/users/erin/RotationProject_Akita/data/genomes/1KG
DIRNAME=$(grep ${INDV} integrated_call_samples_v3.20130502.ALL.panel | awk '{print $3"_"$2"_"$4"_"$1}')
mkdir -p ${DIRNAME}
cd ${DIRNAME}
# mkdir -p ${INDV}
# cd ${INDV}

bash /dors/capra_lab/users/erin/RotationProject_Akita/scripts/makeFASTA/make.genome.1000kg.indiv.sh ${CHR} ${INDV} ${VCFPTH} ${PFX} ${SFX}