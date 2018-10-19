#processing new raw fastq files passed along by Lee Stanish October 2018.
#Raw sequences were joined using fastq-join.
#Seqs demultiplexed using split_libraries_fastq.py, minimal filtering as suggested by dada2 tutorial.
#Details in evernote notebook. Unfortunately, processing these required a lot of "by-hand" steps.
#From here we have per-sample fastq files to be processed in dada2.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
library(dada2)
source('NEFI_functions/tic_toc.r')

#Set path to sequences
path <- 
