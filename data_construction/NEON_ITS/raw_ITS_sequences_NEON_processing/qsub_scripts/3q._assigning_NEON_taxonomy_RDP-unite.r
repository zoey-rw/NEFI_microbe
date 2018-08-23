#!/bin/bash -l
#This script is for Downloading all sequences associated with Bahram et al. 2018 fromt the SRA. 235 samples, 470 files.
#
########################################
####      commands for scc qsub     ####
########################################
#Specfiy hard time limit for the job. 100 hours is plenty.
#$ -l h_rt=100:00:00
#run on a few processors (dada2 RDP classifier can be run in parallel)
#$ -pe omp 8
#$ -l mem_total=510G #huge memory job.
#
#Give the job a name
#$ -N NEON_tax
#
# Merge stderr into the stdout file, to reduce clutter
#$ -j y
#$ -o $JOB_NAME.log
#
#
# Have the system send mail when the job begins and when the job is aborted or ended
#$ -m ae
#
# end of qsub arguments
#
########################################
#### begin commands to run R script ####
########################################
#
#
# load necessary modules 
module load R/3.4.0
module load python/2.7.7
module load qiime/1.9.0
#
#move to directory that contains the project
cd /projectnb/talbot-lab-data/caverill/NEFI_microbe
#
# in the directory specified above, invoke this function:
Rscript data_construction/NEON_ITS/raw_ITS_sequences_NEON_processing/3._assigning_NEON_taxonomy_RDP-unite.r
#
#
#End of commands.
#
