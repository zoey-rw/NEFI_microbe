#!/bin/bash -l
# qsun for fitting functional group models.
########################################
####      commands for scc qsub     ####
########################################
#Specfiy hard time limit for the job.
#$ -l h_rt=120:00:00
#
#Use N processors on a single machine. Running 3 chains in parallel, but i think a 4th will make things move better.
#$ -pe omp 28
#
#Give the job a name
#$ -N bahram_dmulti_JAGS
#
# Merge stderr into the stdout file, to reduce clutter
#$ -j y
#$ -o $JOB_NAME.log
#
# Request buyin nodes
#$ -l buyin
#
# Have the system send mail when the job begins and when the job is aborted or ended
#$ -m ae
#
# Inherit the current environment (load modules python/2.7.7, qiime, and find binaries)
# Make sure th load those modules in the command line before you submit the qsub
#$ -V 
#
# end of qsub arguments
#
########################################
#### begin commands to run R script ####
########################################
#
#
# load necessary modules 
module load R/3.4.3
module load jags
#
# in the directory specified above, invoke this function:
Rscript 16S/spatial_analysis/dmulti-ddirch/1._prior_dmulti.ddirch_fit_16S.r
#
#
#End of commands.
#