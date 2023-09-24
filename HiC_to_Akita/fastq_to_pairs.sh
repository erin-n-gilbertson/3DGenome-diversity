#!/bin/bash
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -l mem_free=2G     # job requires up to 1 GiB of RAM per slot
#$ -l h_rt=72:00:00   # job requires up to 24 hours of runtime

load_conda
conda activate hic
module load CBI fastqc bwa samtools 
#conda activate modern3d

fastq1=HG00864_GCAATTCC-CCACAACA_HF35FDSXY_L002_001.R1.fastq.gz
fastq2=HG00864_GCAATTCC-CCACAACA_HF35FDSXY_L002_001.R2.fastq.gz
prefix=HG00864

outdir=/wynton/group/capra/projects/modern_human_3Dgenome/data/experimental

cd $outdir
mkdir QC_R1
mkdir QC_R2

#fastqc -t ${NSLOTS} $fastq1 -o $outdir/QC_R1
#fastqc -t ${NSLOTS} $fastq2 -o $outdir/QC_R2
#bwa index /wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/hg38_reference.fa

#bwa mem -t ${NSLOTS} -SP5M /wynton/group/capra/projects/modern_human_3Dgenome/data/genomes/hg38_reference.fa $fastq1 $fastq2 | samtools view -Shb - > $prefix.bam
samtools view -h $prefix.bam | {
	pairtools parse -c /wynton/group/capra/projects/modern_human_3Dgenome/data/experimental/hg38.chromsizes --add-columns mapq
} | {
	pairtools sort --nproc ${NSLOTS} --memory 32G --compress-program lz4c --tmpdir $outdir --output $prefix.sam.pairs.gz
}
pairtools dedup --mark-dups --output-dups - --output-unmapped - --output $prefix.marked.sam.pairs.gz $prefix.sam.pairs.gz
pairix $prefix.marked.sam.pairs.gz
pairtools stats --cmd-in 'pbgzip -dc -n '${NSLOTS}'' -o $prefix.marked.pairs.stats $prefix.marked.sam.pairs.gz
pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' --output-rest $prefix.unmapped.sam.pairs.gz --output temp.gz $prefix.marked.sam.pairs.gz
pairtools split --output-pairs temp1.gz temp.gz
pairtools select 'True' --chrom-subset /wynton/group/capra/projects/modern_human_3Dgenome/data/experimental/hg38.chromsizes -o $prefix.dedup.pairs.gz temp1.gz
pairix $prefix.dedup.pairs.gz

