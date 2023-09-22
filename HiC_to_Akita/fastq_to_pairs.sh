#!/bin/bash


module load CBI fastqc bwa samtools 
#conda activate modern3d

fastq1=$1
fastq2=$2
prefix=$3
nThreads=$4
outdir=$5

cd $outdir
mkdir QC_R1
mkdir QC_R2

#fastqc -t $nThreads $fastq1 -o $outdir/QC_R1
#fastqc -t $nThreads $fastq2 -o $outdir/QC_R2
#bwa index /wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/hg38_reference.fa

bwa mem -t $nThreads -SP5M /wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/hg38_reference.fa $fastq1 $fastq2 | samtools view -Shb - > $prefix.bam
samtools view -h $prefix.bam | {
	pairtools parse -c /wynton/group/capra/projects/modern_human_3Dgenome/data/experimental/hg38.chromsizes --add-columns mapq
} | {
	pairtools sort --nproc $nThreads --memory 32G --compress-program lz4c --tmpdir $outdir --output $prefix.sam.pairs.gz
}
pairtools dedup --mark-dups --output-dups - --output-unmapped - --output $prefix.marked.sam.pairs.gz $prefix.sam.pairs.gz
pairix $prefix.marked.sam.pairs.gz
pairtools stats --cmd-in 'pbgzip -dc -n '$nThreads'' -o $prefix.marked.pairs.stats $prefix.marked.sam.pairs.gz
pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' --output-rest $prefix.unmapped.sam.pairs.gz --output temp.gz $prefix.marked.sam.pairs.gz
pairtools split --output-pairs temp1.gz temp.gz
pairtools select 'True' --chrom-subset /wynton/group/capra/projects/modern_human_3Dgenome/data/experimental/hg38.chromsizes -o $prefix.dedup.pairs.gz temp1.gz
pairix $prefix.dedup.pairs.gz

