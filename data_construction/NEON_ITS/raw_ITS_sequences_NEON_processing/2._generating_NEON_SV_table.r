#Processing raw ITS sequences from NEON.
#1. These are forward reads only, reverse was discarded.
#2. Quality filtering has already been done.
#3. Go through and separate reads by sample, contruct per sample files.
#4. Make SV tables from per sample files.
#5. remove chimeras, save SV table.

#clear environment, load functions and packages.
rm(list=ls())
source('paths.r')
library(data.table)

#get file paths
files <- list.files(NEON_ITS.dir)
files <- files[grep('.fna',files)]
