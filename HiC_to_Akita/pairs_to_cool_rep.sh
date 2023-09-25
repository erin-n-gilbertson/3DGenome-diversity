
#!/bin/bash
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -l mem_free=2G     # job requires up to 1 GiB of RAM per slot
#$ -l h_rt=72:00:00   # job requires up to 24 hours of runtime


module load CBI miniconda3/4.12.0-py39
conda activate hic
module load CBI fastqc bwa samtools 


prefix=HG00864

outdir=/wynton/group/capra/projects/modern_human_3Dgenome/data/experimental

cd $outdir

cooler cload pairix -s 2 --assembly hg38 -p ${NSLOTS} hg38.chromsizes:2048 $prefix.dedup.pairs.gz $prefix.cool
cooler balance --max-iters 500 --convergence-policy store_final $prefix.cool 
#cooler zoomify -n ${NSLOTS} --balance --balance-args '--max-iters 500 --convergence-policy store_final' -r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000 -o $prefix.mcool $prefix.cool
