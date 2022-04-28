#!/bin/bash
#$ -m a
#$ -o /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/launch.o
#$ -e /wynton/group/capra/projects/modern_human_3Dgenome/stdout/make_genome/launch.e


qsub /wynton/group/capra/projects/modern_human_3Dgenome/bin/makeFASTA/launch.make.genome.1KG_afr.q.sh
# qsub /wynton/group/capra/projects/modern_human_3Dgenome/bin/makeFASTA/launch.make.genome.1KG_amr.q.sh
# qsub /wynton/group/capra/projects/modern_human_3Dgenome/bin/makeFASTA/launch.make.genome.1KG_eas.q.sh
# qsub /wynton/group/capra/projects/modern_human_3Dgenome/bin/makeFASTA/launch.make.genome.1KG_eur.q.sh
# qsub /wynton/group/capra/projects/modern_human_3Dgenome/bin/makeFASTA/launch.make.genome.1KG_sas.q.sh


#

#
