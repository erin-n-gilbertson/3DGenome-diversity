This directory contains a scripts to generate 3d maps from FASTAs using Akita. 
- config files can be customized for file paths and lists of FASTA input files as needed
- runAkita.sh: SGE job script that parses config file and submits an task array per individual. Should add '-t' SGE parameter to header to indicate number of jobs in array job.
- runAkita.one_individual.py: launched by above, python script to make Akita predictions for one individual given a list of genome chunks which can be specified in the config file. 
