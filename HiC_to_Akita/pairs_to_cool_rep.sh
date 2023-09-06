#!/bin/bash

prefix=$1
nThreads=$2
outdir=$3

cd $outdir
cooler cload pairix -s 2 --assembly hg38 -p $nThreads /wynton/home/pollard/skuang/4DN/HiC/Data/GRCh38_EBV.chrom.sizes:1000 $prefix.dedup.pairs.gz $prefix.cool
cooler balance $prefix.cool
cooler zoomify -n $nThreads --balance --balance-args '--max-iters 500 --convergence-policy store_final' -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 -o $prefix.mcool $prefix.cool
