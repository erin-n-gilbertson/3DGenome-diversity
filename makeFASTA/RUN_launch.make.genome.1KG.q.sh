#!/bin/bash
#$ -m a
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/launch.o
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/launch.e


qsub launch.make.genome.1KG_afr.q.sh
qsub launch.make.genome.1KG_amr.q.sh
qsub launch.make.genome.1KG_eas.q.sh
qsub launch.make.genome.1KG_eur.q.sh
qsub launch.make.genome.1KG_sas.q.sh


#

#
