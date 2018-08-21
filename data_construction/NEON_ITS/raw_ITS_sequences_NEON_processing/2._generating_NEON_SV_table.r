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

#get .fna file paths
seq.path <- NEON_ITS.dir
files <- list.files(seq.path)
files <- files[grep('.fasta',files)]

#Get reverse complement of reverse primers.
#These are ITS1 reads with reverse read discarded. Some sequences have reverse primer, some dont. 
#All have forward primer already trimmed.
rc.rev.primer <- 'GCATCGATGAAGAACGCAGC'


#use bbduk to trim reverse primer if present.
files <- files[1:2] #subset to test
bbduk.path <- 'NEFI_functions/bbmap/bbduk.sh' #path to bbduk function within the bbmap directory.
for(i in 1:length(files)){
  #quality filter fastq files using qiime.
   output.dir <- 'q.trim/'
  sample.name <- files[i]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
   input.path <- paste0(seq.path,sample.name,'.fasta')
  output.path <- paste0(seq.path,output.dir,sample.name,'.fna')
  cmd <- paste0(bbduk.path,
                ' literal=',rc.rev.primer,
                ' ktrim=r k=10',
                ' in=',input.path,
                ' out=',output.path)
  system(cmd)
}

#above code made sequence take up more than one line, which interferes with code I use to make SV table.
#we fix this here.
q.final.dir <- paste0(seq.path,'q.final/')
system(paste0('mkdir -p ',q.final.dir))
for(i in 1:length(files)){
  sample.name <- files[i]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
   input.path <- paste0(seq.path,'q.trim/',sample.name,'.fna')
  output.path <- paste0(       q.final.dir,sample.name,'.fna')
  system(command = paste0("perl -pe '/^>/ ? print \"\n\" : chomp' ",
                          input.path,
                          " | tail -n +2 > ",
                          output.path), intern = TRUE)
}

#remove q.trim
cmd <- paste0('rm -rf ',seq.path,'q.trim')
system(cmd)
