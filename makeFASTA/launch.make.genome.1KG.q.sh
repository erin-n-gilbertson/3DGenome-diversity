#!/bin/bash
#$ -M erin.gilbertson@ucsf.edu
#$ -m abe
#$ -o /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/stdout
#$ -e /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/stdout

## the called SLURM script usees an internal array (ARRY) of chromosome names that the slurm array IDs are then used to index
## ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)



###### TODO_Erin
SAMPLES=(NA19794)


for i in "${SAMPLES[@]}"
do
	INDIV=${i}
	qsub -N make.genome.$INDIV -l mem_free=19G -t 22 -l h_rt=2:00:00 -V /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/bin/makeFASTA/make.genome.1KG.indiv.q
	#qsub -N make.genome.$INDIV -l mem_free=40G -t 1-4 -l h_rt=2:00:00 -V /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/bin/makeFASTA/make.genome.1KG.indiv.q
	#qsub -N make.genome.$INDIV -l mem_free=19G -t 5-23 -l h_rt=2:00:00 -V /wynton/home/capra/egilbertson/projects/modern_human_3Dgenome/bin/makeFASTA/make.genome.1KG.indiv.q

done

#
