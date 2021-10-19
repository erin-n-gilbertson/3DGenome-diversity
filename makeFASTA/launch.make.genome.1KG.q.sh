#!/bin/bash

## the called SLURM script usees an internal array (ARRY) of chromosome names that the slurm array IDs are then used to index
## ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)



###### TODO_Erin
SAMPLES=(NA19794)


for i in "${SAMPLES[@]}"
do
	INDIV=${i}
	qsub --N make.genome.$INDIV -l mem_free=40G -t 0-3 -l h_rt=2:00:00 -V make.genome.1KG.indiv.q
	qsub --N make.genome.$INDIV -l mem_free=19G -t 4-22 -l h_rt=2:00:00 -V make.genome.1KG.indiv.q

done

#
