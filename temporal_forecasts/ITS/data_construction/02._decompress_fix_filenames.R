# 1. Decompress .tar.gz files
# 2. Match .fastq file names to the sequences that they correspond with
# 3. Remove .fastq files that are not from our 5 sites

rm(list=ls())
library(RCurl)
library(neonUtilities)
library(utils)
source('NEFI_functions/tic_toc.r')
source('paths.r')

# set sequence and output directory
seq.dir <- paste0(big_data_dir, 'NEON_2014-2017/ITS/')
#download.dir <- paste0(seq.dir, 'raw_fastq_ITS/')
download.dir <- "/projectnb/talbot-lab-data/NEFI_data/big_data/NEON_2014-2017/ITS/test/" #testing
#output.dir <- paste0(seq.dir, '2_per_sample_fastq.gz/')
output.dir <- "/usr3/graduate/zrwerbin/hpc/home/minardsmitha/NEON/16S_ITS_2017/Sept_12_Run_BFDG8/RAW_FASTQ/ITS/RAW_Upload_to_BOX/R1/" #testing

# go to directory of tar files and decompress all of them
tarfiles <- list.files(download.dir)
system(paste0("cd ", download.dir, '; for f in *.tar.gz; do tar -zxvf "$f" -C ', output.dir, '; done'))

## UNNEST FILES ##
# ...figure out how to do this...

# then, fix file names 
files <- list.files(output.dir)
files <- files[grep('.fastq',files)]
for(i in 1:length(files)){
  this.file <- paste0(output.dir, files[i])
  base.file <- basename(this.file)
  R1.or.R2 <- substr(base.file, nchar(base.file) - 8, nchar(base.file)) 
  df.file.name <- gsub(x = base.file, pattern = "_R1.fastq|_R2.fastq", replacement = "") #remove R1/R2 suffix
  metadata_filenames <- gsub(x = its_metadata$processedSeqFileName, pattern = "_fastq.gz|_fastq|.fastq|.fasta|_fastq.tar.gz", replacement = "") #remove R1/R2 suffix
  metadata_dnaSampleID <- its_metadata[which(metadata_filenames == df.file.name),]$dnaSampleID # match within metadata
  if (length(metadata_dnaSampleID)==0) next # if there's no match, move on
  new.name <- paste0(metadata_dnaSampleID, R1.or.R2)
  cmd <- paste("mv", this.file, paste0(output.dir, new.name))
  print(cmd) #testing
  #system(cmd)
}
