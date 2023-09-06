#!/bin/bash

fastq1=$1
fastq2=$2
prefix=$3
nThreads=$4
outdir=$5

cd $outdir
mkdir QC_R1
mkdir QC_R2

fastqc -t $nThreads $fastq1 -o $outdir/QC_R1
fastqc -t $nThreads $fastq2 -o $outdir/QC_R2
bwa mem -t $nThreads -SP5M /wynton/group/gladstone/references/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta $fastq1 $fastq2 | samtools view -Shb - > $prefix.bam
samtools view -h $prefix.bam | {
	pairtools parse -c /wynton/home/pollard/skuang/4DN/HiC/Data/GRCh38_EBV.chrom.sizes --add-columns mapq
} | {
	pairtools sort --nproc $nThreads --memory 32G --compress-program lz4c --tmpdir $outdir --output $prefix.sam.pairs.gz
}
pairtools dedup --mark-dups --output-dups - --output-unmapped - --output $prefix.marked.sam.pairs.gz $prefix.sam.pairs.gz
pairix $prefix.marked.sam.pairs.gz
pairtools stats --cmd-in 'pbgzip -dc -n '$nThreads'' -o $prefix.marked.pairs.stats $prefix.marked.sam.pairs.gz
pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' --output-rest $prefix.unmapped.sam.pairs.gz --output temp.gz $prefix.marked.sam.pairs.gz
pairtools split --output-pairs temp1.gz temp.gz
pairtools select 'True' --chrom-subset /wynton/home/pollard/skuang/4DN/HiC/Data/GRCh38_EBV.chrom.sizes -o $prefix.dedup.pairs.gz temp1.gz
pairix $prefix.dedup.pairs.gz

