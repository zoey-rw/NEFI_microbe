#Constructing SV table froM Tedersoo 2014 .fastq files downloaded from the SRA.
#this is using qiime quality filtering, then trimmomatic to remove primers.
#After this we de-replicate sequences, construct SV table, remove chimeras using dada2.
#assign taxonomy using RDP via dada2.
#starting with two test files.
#require qsub script to load the following modules: R/3.4.0, python/2.7.7, qiime/1.9.0
#clear environment, load packages.
rm(list=ls())
source('paths.r')
library(data.table)

#begin by testing this script with 2 .fastq files from the tedersoo study.
seq.path <- ted.seq.dir
#seq.path <- '/projectnb/talbot-lab-data/caverill/ted_test_fastq/' #for testing.
#seq.path <- '/fs/data3/caverill/NEFI_microbial_data/ted_test_fastq/'

#reverse primers (there is a flex position)
rev.primers <- 'TCCTGCGCTTATTGATATGC,TCCTCCGCTTATTGATATGC'
#foward primers: there are 6. These are their reverse complements.
rc.fwd.primers <- ('CAGCGTTCTTCATCGATGACGAGTCTAG,CTGCGTTCTTCATCGTTGACGAGTCTAG,CTGCGTTCTTCATCGGTGACGAGTCTAG,CTACGTTCTTCATCGATGACGAGTCTAG,CCACGTTCTTCATCGATGACGAGTCTAG,CAGCGTTCTTCATCGATGACGAGTCTAG')

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

#update fastq list in case a file or two didn't make ith through split_libraries_fastq.py.
fastq.files <- list.files(paste0(seq.path,'q.filter/'))
fastq.files <- fastq.files[grep('.fna',fastq.files)]
fastq.files <- gsub('.fna','.fastq',fastq.files)

#from here we need to trim out primers and leading adapter/barcode sequence which is variable length.
#DOE has a great tool to find a primer and trim anything preceding in, called bbuk.sh in the bbmap package.
#Find theat here: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/
#I copied this indivudal script into the NEFI_tools directory.
#This also trim 3' "forward" primer, and requires sequences to be 100bp long.
for(i in 1:length(fastq.files)){
  #quality filter fastq files using qiime.
  sample.name <- fastq.files[i]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
  output.dir1 <- 'q.trim.L/'
  bbduk.path <- 'NEFI_functions/bbmap/bbduk.sh'
  cmd <- paste0(bbduk.path,
                ' literal=',rev.primers,
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
cmd <- paste0('rm -rf ',seq.path,'q.filter')
system(cmd)

#above code made sequence take up more than one line, which interferes with code I use to make SV table.
#we fix this here.
system(paste0('mkdir -p ',seq.path,'q.final'))
for(i in 1:length(fastq.files)){
  sample.name <- fastq.files[i]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
  system(command = paste0("perl -pe '/^>/ ? print \"\n\" : chomp' ",
                          seq.path,"q.trim/",sample.name,
                          ".fna | tail -n +2 > ",
                          seq.path,"q.final/",sample.name,".fna"), intern = TRUE)
}

#de-replicate and generate SV table.
#### Build an ASV table from all sequence files. ####
cat(paste0('Building ASV tables...\n'))
for(j in 1:length(fastq.files)){
  sample.name <- fastq.files[j]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
  #get just the sequences into a separate file, without additional shit.
  
    file <- paste0(seq.path,'q.final/'    ,sample.name,'.fna')
  s.file <- paste0(seq.path,'q.final/seq.',sample.name,'.fna')
  cmd <- paste0("sed -n '0~2p' ",file,' > ',s.file)
  system(cmd)
  
  #Convert sequences to a sequence table. Write table to a file. No need to load all sequences into R memory.
  c.file <- paste0(seq.path,'counts.',sample.name,'.txt')   #sample-specific ASV table.
  pre <- paste0('cat ',s.file,' | sort | uniq -c > ', c.file)
  system(pre)
  asv <- data.table::data.table(read.table(c.file))
  colnames(asv) <- c(sample.name,'seq')
  data.table::setkey(asv,seq)
  
  #merge multiple sample-specific ASV tables.
  if(j == 1){ out = asv}
  if(j > 1){ out <- merge(out, asv, all = T)} #this merge really requires data.table be loaded into the environment.
  
  #remove duplicated seq and counts files.
  system(paste0('rm ',c.file))
  system(paste0('rm ',s.file))
}

#convert back to dataframe, replace NA values with zeros.
out <- as.data.frame(out)
out[is.na(out)] <- 0

#transpose to be consistent with dada2
#this is being weird.
t.out <- t(out[,-1])
colnames(t.out) <- as.character(out[,1])

#convert from numeric dataframe to integer matrix. this is important for dada2 commands downstream.
t.out <- as.matrix(t.out)
t.out <- apply (t.out, c (1, 2), function (x) {(as.integer(x))})
cat('ASV table built!\n')

#Remove chimeras.
cat('Removing chimeras...\n')
t.out_nochim <- dada2::removeBimeraDenovo(t.out, method = 'consensus', multithread = T)
cat('Chimeras removed.\n')

#save output.
output_filepath <- paste0(seq.path,'SV_table.rds')
saveRDS(t.out_nochim, output_filepath)
cat('script complete.\n')


#forward primers 5' -> 3' (before reverse complement), which are actually at the ends of these reads.
#CTAGACTCGTCATCGATGAAGAACGCAG
#CTAGACTCGTCAACGATGAAGAACGCAG
#CTAGACTCGTCACCGATGAAGAACGCAG
#CTAGACTCGTCATCGATGAAGAACGTAG
#CTAGACTCGTCATCGATGAAGAACGTGG
#CTAGACTCGTCATCGATGAAGAACGCTG
