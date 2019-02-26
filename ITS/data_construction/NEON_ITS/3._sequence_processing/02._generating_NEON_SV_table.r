#Processing raw ITS sequences from NEON. This script takes ~3 hours to run with 36 cores in parallel, big memory.
#Slowest part is the merging on the sample-specific SV tables.
#Time is taken scales linearly in log-log space with number of samples merged.
#Someone clever could write a hierarchical merging algorithm that might speed this way up, resulting in fewer merging comparisons.

#1. clear environment, load functions and packages.----
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(data.table)
library(doParallel)

#register parallel environment.
n <- detectCores()
registerDoParallel(cores=n)

#2. Setup paths.----
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


#3. use bbduk to trim reverse primer if present.----
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

#4. collapse to normal format.----
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

#5. de-replicate and generate SV table.----
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

#Build an ASV table from all sequence files.
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

#6. Merge sample specific SV tables.----
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
  if(to_check == round(to_check)){cat(j,'of',length(files),'samples merged. ');toc()}
  
}
cat('Sample-specific SV tables merged. ')
toc()

#7. dial in the SV table to work with dada2.----
cat('Dialing in merged SV table...\n')
#replace NA values with zeros.
#You need these data.table trick when the SV table is BIG.
DT.for.set.sqln  <- function(x) { for (j in seq_len(ncol(x)))
  set(x,which(is.na(x[[j]])),j,0) }
out <- data.table(out)
DT.for.set.sqln(out)

#transpose to be consistent with dada2. These steps are fine.
t.out <- t(out[,-1])
colnames(t.out) <- as.character(out[,1])

#convert from numeric dataframe to integer matrix. this is important for dada2 commands downstream.
t.out <- as.matrix(t.out)
mode(t.out) <- "integer"

#drop singletons and reads <100bp.
t.out <- t.out[,!colSums(t.out) == 1]
#remove reads less than 100bp
t.out <- t.out[,nchar(colnames(t.out)) > 99]
cat('ASV table dialed!\n')
SV_pre.chimera.path <- paste0(seq.path,'SV_pre.chimera_table.rds')
saveRDS(t.out,SV_pre.chimera.path)


#clean up (delete) q.final sample-specific sv. directories.
system(paste0('rm -rf ',seq.path,'q.final'))
system(paste0('rm -rf ',sv.dir.path))

#8. Remove chimeras using dada2.----
######## We may be able to return this here now that I filter out singletons and reads < 100bp
#cat('Removing chimeras...\n')
#tic()
#t.out_nochim <- dada2::removeBimeraDenovo(t.out, method = 'consensus', multithread = T)
#cat('Chimeras removed.\n')
#toc()

#9. Final save and cleanup.----
#save output.
#output_filepath <- paste0(seq.path,'SV_table.rds')
#saveRDS(t.out_nochim, output_filepath1)
#saveRDS(t.out_nochim, output_filepath2)

cat('script complete.\n')
