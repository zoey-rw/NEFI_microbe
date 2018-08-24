#Processing raw ITS sequences from NEON.
#1. These are forward reads only, reverse was discarded.
#2. Quality filtering has already been done.
#3. Go through and separate reads by sample, contruct per sample files.
#4. Make SV tables from per sample files.
#5. remove chimeras, save SV table.

#clear environment, load functions and packages.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(data.table)
library(doParallel)

#register parallel environment.
n <- detectCores()
registerDoParallel(cores=n)

#get sequence file paths
seq.path <- NEON_ITS.dir
files <- list.files(seq.path)
files <- files[grep('.fasta',files)]

#output file path.
output_filepath1 <-  paste0(seq.path,'SV_table.rds')
output_filepath2 <- NEON_SV.table.path

#Get reverse complement of reverse primers.
#These are ITS1 reads with reverse read discarded. Some sequences have reverse primer, some dont. 
#All have forward primer already trimmed.
rc.rev.primer <- 'GCATCGATGAAGAACGCAGC'


#use bbduk to trim reverse primer if present.
#files <- files[1:2] #subset to test
cat('Trimming primers...\n')
tic()
bbduk.path <- 'NEFI_functions/bbmap/bbduk.sh' #path to bbduk function within the bbmap directory.
#for(i in 1:length(files)){
foreach(i = 1:length(files)) %dopar% {
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
#for(i in 1:length(files)){
foreach(i = 1:length(files)) %dopar% {
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
cat('Primers trimmed!\n')
toc()

#de-replicate and generate SV table.
#Some files have zero lines. remove these.
path_to_check <- paste0(seq.path,'q.final/')
cmd <- paste0('find ',path_to_check,' -size 0 -delete')
system(cmd)
#update file list.
files <- list.files(paste0(seq.path,'q.final/'))
files <- files[grep('.fna',files)]
files <- gsub('.fna','.fastq',files)

#make sample-specific SV field irectiory.
sv.dir.path <- paste0(seq.path,'sv.files/')
cmd <- paste0('mkdir -p ',sv.dir.path)
system(cmd)

#### Build an ASV table from all sequence files. ####
cat(paste0('Building sample specific SV tables...\n'))
tic()
#for(j in 1:length(files)){
foreach(j = 1:length(files)) %dopar% {
  sample.name <- files[j]
  sample.name <- substr(sample.name,1,nchar(sample.name)-6)
  #get just the sequences into a separate file, without header lines.
  file <- paste0(seq.path,'q.final/'      ,sample.name,'.fna')
  s.file <- paste0(seq.path,'q.final/seq.',sample.name,'.fna')
  cmd <- paste0("sed -n '0~2p' ",file,' > ',s.file)
  system(cmd)
  
  #Convert sequences to a sequence table. Write table to a file. No need to load all sequences into R memory.
  c.file <- paste0(sv.dir.path,'sv.',sample.name,'.txt')   #sample-specific ASV table.
  pre <- paste0('cat ',s.file,' | sort | uniq -c > ', c.file)
  system(pre)
  
  #clean up s.file
  system(paste0('rm ',s.file))
}
cat('Sample-specific SV tables constructed.\n')
toc()

#update SV files.
files <- list.files(sv.dir.path)

#actually merge the files one at a time. This is the slow part, the big memory part.
cat('Merging sample-specific SV files...\n')
tic()
for(j in 1:length(files)){
  sample.name <- files[j]
  sample.name <- substr(sample.name,1,nchar(sample.name)-4)
  c.file <- paste0(sv.dir.path,files[j])
  asv <- data.table(read.table(c.file))
  colnames(asv) <- c(sample.name,'seq')
  setkey(asv,seq)
  
  #merge multiple sample-specific SV tables.
  if(j == 1){ out = asv}
  if(j  > 1){ out <- merge(out, asv, all = T)} #this merge really requires data.table be loaded into the environment.
  
  #report every 100 samples merged and upadte how long its been running.
  to_check <- j/100
  if(to_check == round(to_check)){cat(j,'of',length(files),'samples merged.\n')}
  toc()
}
cat('Sample-specific SV tables merged. total time:')
toc()

cat('Dialing in merged SV table...\n')
#convert back to dataframe, replace NA values with zeros.
out <- as.data.frame(out)
out[is.na(out)] <- 0

#transpose to be consistent with dada2.
t.out <- t(out[,-1])
colnames(t.out) <- as.character(out[,1])

#convert from numeric dataframe to integer matrix. this is important for dada2 commands downstream.
t.out <- as.matrix(t.out)
t.out <- apply (t.out, c (1, 2), function (x) {(as.integer(x))})
cat('ASV table dialed!\n')

#clean up (delete) q.final sample-specific sv. directories.
#system(paste0('rm -rf ',seq.path,'q.final'))
#system(paste0('rm -rf ',sv.dir.path))

#Remove chimeras.
cat('Removing chimeras...\n')
tic()
t.out_nochim <- dada2::removeBimeraDenovo(t.out, method = 'consensus', multithread = T)
cat('Chimeras removed.\n')
toc()

#sequences must be at least 100bp.
t.out_nochim <- t.out_nochim[,nchar(colnames(t.out_nochim)) > 99]

#save output.
output_filepath <- paste0(seq.path,'SV_table.rds')
saveRDS(t.out_nochim, output_filepath1)
saveRDS(t.out_nochim, output_filepath2)
cat('script complete.\n')
