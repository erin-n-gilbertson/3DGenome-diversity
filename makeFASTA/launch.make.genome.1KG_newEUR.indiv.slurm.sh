#!/bin/bash
#SBATCH --reservation=AtRiskUnstableEnvironment
## the called SLURM script usees an internal array (ARRY) of chromosome names that the slurm array IDs are then used to index
## ARRY=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

########################
### Test run:

# INDIV='NA19159'	# AFR YRI female
# VCFPATH='/gpfs51/dors2/capra_lab/data/1000_genomes_project/phase3/vcf/'
# PREFIX='ALL.chr'
# SUFFIX='.phase3_shapeit2_mvncall_integrated_*.genotypes.vcf.gz'

# sbatch --job-name=make.genome.$INDIV --mem-per-cpu=40000 --array=0-3 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm
# sbatch --job-name=make.genome.$INDIV --mem-per-cpu=19000 --array=4-22 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm

#########
## Erin Test

# INDIV='NA18624'	
# VCFPATH='/gpfs51/dors2/capra_lab/data/1000_genomes_project/phase3/vcf/'
# PREFIX='ALL.chr'
# SUFFIX='.phase3_shapeit2_mvncall_integrated_*.genotypes.vcf.gz'

# sbatch --job-name=make.genome.$INDIV --mem-per-cpu=40000 --array=0-3 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm
# sbatch --job-name=make.genome.$INDIV --mem-per-cpu=19000 --array=4-22 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm

#######################
### ROUND 1 Select random individuals from populations:

# [rinkerd@capra1 1KG.individuals]$ grep AFR integrated_call_samples_v3.20130502.ALL.panel | grep female | shuf -n 5
# HG03212	MSL	AFR	female
# NA18870	YRI	AFR	female
# NA19378	LWK	AFR	female
# HG03539	GWD	AFR	female
# HG03105	ESN	AFR	female
# [rinkerd@capra1 1KG.individuals]$ grep EAS integrated_call_samples_v3.20130502.ALL.panel | grep female | shuf -n 5
# HG00560	CHS	EAS	female
# NA18595	CHB	EAS	female
# HG01851	KHV	EAS	female
# NA19077	JPT	EAS	female
# HG00978	CDX	EAS	female
# [rinkerd@capra1 1KG.individuals]$ grep EUR integrated_call_samples_v3.20130502.ALL.panel | grep female | shuf -n 5
# NA12006	CEU	EUR	female
# HG00285	FIN	EUR	female
# NA20795	TSI	EUR	female
# HG01519	IBS	EUR	female
# HG00261	GBR	EUR	female

###### TODO_Erin
SAMPLES=(HG01773 HG01516 HG01632 HG01668 NA12717 NA12815 NA12156 NA12749 HG00288 HG00323 HG00304 HG00344 NA20821 HG00231 NA20508 NA20797 NA20771 HG00106 HG00122 HG00099)

for i in "${SAMPLES[@]}"
do
	INDIV=${i}
	VCFPATH='/gpfs51/dors2/capra_lab/data/1000_genomes_project/phase3/vcf/'
	PREFIX='ALL.chr'
	SUFFIX='.phase3_shapeit2_mvncall_integrated_*.genotypes.vcf.gz'

	sbatch --job-name=make.genome.$INDIV --mem-per-cpu=40000 --array=0-3 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm
	sbatch --job-name=make.genome.$INDIV --mem-per-cpu=19000 --array=4-22 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm

done

### ROUND1 Clean up ACCRE's problem jobs:
# VCFPATH='/gpfs51/dors2/capra_lab/data/1000_genomes_project/phase3/vcf/'
# PREFIX='ALL.chr'
# SUFFIX='.phase3_shapeit2_mvncall_integrated_*.genotypes.vcf.gz'

# INDIV=HG00978
# sbatch --job-name=make.genome.$INDIV --exclude=cn458,cn1384 --mem-per-cpu=19000 --array=4-6 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm

# INDIV=NA19077
# sbatch --job-name=make.genome.$INDIV --exclude=cn458,cn1384 --mem-per-cpu=19000 --array=20-22 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm

# INDIV=HG03212
# sbatch --job-name=make.genome.$INDIV --exclude=cn458,cn1384 --mem-per-cpu=19000 --array=21 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm

# INDIV=NA12006
# sbatch --job-name=make.genome.$INDIV --exclude=cn458,cn1384 --mem-per-cpu=19000 --array=11 --time=2:00:00 --export=INDV=$INDIV,VCFPTH=$VCFPATH,PFX=$PREFIX,SFX=$SUFFIX make.genome.1KG.indiv.slurm

#######################
### ROUND 2 What Next???