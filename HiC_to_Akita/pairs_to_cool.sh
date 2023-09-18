#!/bin/bash

rep1=$1
rep2=$2
prefix=$3
nThreads=$4
outdir=$5

cd $outdir
ln -s /wynton/group/capra/projects/modern_human_3Dgenome/bin/HiC_to_Akita/run_merge_pairs.sh 
sh run_merge_pairs.sh $prefix.comb $rep1 $rep2
cooler cload pairix -s 2 --assembly hg38 -p $nThreads /wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/hg38.chrom.bed:1000 $prefix.comb.pairs.gz $prefix.cool
cooler balance $prefix.cool
cooler zoomify -n $nThreads --balance --balance-args '--max-iters 500 --convergence-policy store_final' -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 -o $prefix.mcool $prefix.cool


