#!/bin/bash

pairfile=$1
outdir=$2
scriptdir=$3
samplename=$4

python $scriptdir/pairsqc.py -p $pairfile -c /wynton/home/pollard/skuang/4DN/HiC/Data/GRCh38_noM.chrom.sizes -tP -s $samplename -O $samplename -M 8.4
Rscript $scriptdir/plot.r 4 $outdir
