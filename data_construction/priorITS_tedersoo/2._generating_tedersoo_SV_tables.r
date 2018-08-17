#Constructing SV table froM Tedersoo 2014 .fastq files downloaded from the SRA.
#this is using qiime quality filtering, then trimmomatic to remove primers.
#After this we de-replicate sequences, construct SV table, remove chimeras using dada2.
#assign taxonomy using RDP via dada2.
#starting with two test files.
#require qsub script to load the following modules: R/3.4.0, python/2.7.7, qiime/1.9.0
#clear environment, load packages.
rm(list=ls())
source('paths.r')

#begin by testing this script with 2 .fastq files from the tedersoo study.
seq.path <- ted.seq.dir
seq.path <- '/projectnb/talbot-lab-data/caverill/ted_test_fastq/' #for testing.

#reverse primers (there is a flex position)
rev.primers <- 'TCCTGCGCTTATTGATATGC,TCCTCCGCTTATTGATATGC'
#foward primers: there are 6. These are their reverse complements.
rc.fwd.primers <- 'CAGCGTTCTTCATCGATGACGAGTCTAG,CTGCGTTCTTCATCGTTGACGAGTCTAG,CTGCGTTCTTCATCGGTGACGAGTCTAG,CTACGTTCTTCATCGATGACGAGTCTAG,CCACGTTCTTCATCGATGACGAGTCTAG,CAGCGTTCTTCATCGATGACGAGTCTAG'

#get fastq file names, only include files that end in .fastq.
fastq.files <- list.files(seq.path)
fastq.files <- fastq.files[grep('.fastq',fastq.files)]

#loop over all files, perform quality filter using qiime.
for(i in 1:length(fastq.files)){
  #quality filter fastq files using qiime.
  sample.name <- fastq.files[i]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
  output.dir <- 'q.filter/'
  cmd <- paste0('split_libraries_fastq.py -i ',seq.path,sample.name,'.fastq -o ',seq.path,output.dir,
                ' --barcode_type not-barcoded --phred_offset 33 --sample_ids ',sample.name)
  system(cmd)
  #rename seqs.fna file.
  cmd <- paste0('mv ',seq.path,output.dir,'seqs.fna ',seq.path,output.dir,sample.name,'.fna')
  system(cmd)
  #remove log and txt files.
  cmd <- paste0('rm ',seq.path,output.dir,'*.txt')
  system(cmd)
}

#from here we need to trim out primers and leading adapter/barcode sequence which is variable length.
#DOE has a great tool to find a primer and trim anything preceding in, called bbuk.sh in the bbmap package.
#Find theat here: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
#This also trim 3' "forward" primer, and requires sequences to be 100bp long.
for(i in 1:length(fastq.files)){
  #quality filter fastq files using qiime.
  sample.name <- fastq.files[i]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
  output.dir1 <- 'q.trim.L/'
  cmd <- paste0('/projectnb/talbot-lab-data/caverill/bbmap/bbduk.sh ',
                'literal=',rev.primers,
                ' ktrim=l k=10 minlen=100 ',
                'in=',seq.path,'q.filter/',sample.name,'.fna out=',seq.path,output.dir1,sample.name,'.fna')
  system(cmd)
  #now trim 3' end since all "forward" primers are 28bp long.
  output.dir2 <- 'q.trim.R/'
  cmd <- paste0('/projectnb/talbot-lab-data/caverill/bbmap/bbduk.sh ',
                'literal=',rc.fwd.primers,
                ' ktrim=r k=10 ',
                'in=',seq.path,output.dir1,sample.name,'.fna out=',seq.path,output.dir2,sample.name,'.fna')
  system(cmd)
}

#clean up and rename some things.
cmd <- paste0('rm -rf ',seq.path,'q.trim.L')
system(cmd)
cmd <- paste0('mv ',seq.path,'q.trim.R ',seq.path,'q.trim')
system(cmd)
#cmd <- 'rm -rf ',seq.path,'q.filter' #once you are sure everything is working do this.
#system(cmd)

#de-replicate and generate SV table.







#forward primers 5' -> 3' (before reverse complement), which are actually at the ends of these reads.
#CTAGACTCGTCATCGATGAAGAACGCAG
#CTAGACTCGTCAACGATGAAGAACGCAG
#CTAGACTCGTCACCGATGAAGAACGCAG
#CTAGACTCGTCATCGATGAAGAACGTAG
#CTAGACTCGTCATCGATGAAGAACGTGG
#CTAGACTCGTCATCGATGAAGAACGCTG
