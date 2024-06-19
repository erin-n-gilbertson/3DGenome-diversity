This directory contains the scripts and config file to generate haploid FASTA sequences for input into Akita. Scripts are listed here in the order they should run. 
Input is VCFs which can be downloaded from the 1KG ftp site: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/

- config_makeFASTA_1KG_SNVs.ini: configuration file with all of the necessary parameters for running the following scripts. Any changes to be made should be made here. File paths can be specified to reflect the users project structure. 
- launch.make.genome.1KG.sh: bash job script to be run as an array job for the number of individuals. For each individual launches array jobs over chromosomes. 
- make.genome.1KG.one_indiv.sh: Launched by above;  bash job script that launches bash scripts to create FASTA files per chromosome.
- make.genome.1KG.one_chromosome.ch: bash script that does the bulk of the work in this process. For a given chromosome for a given individual creates haploid FASTA file keeping every alternate allele. 
