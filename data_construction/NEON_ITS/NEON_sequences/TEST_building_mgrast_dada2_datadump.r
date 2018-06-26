#downloading a NEON site-month set of sequences, processing in DADA2 to ASV table and RDP taxonomy table.
#chimera checking takes a super lnog time. Make a separate function?
#RDP taxonomy should probably be a separate function.
rm(list=ls())

#ITS ASV master path. Function should sepcify ITS or 16S folder.
main.path <- '/fs/data3/caverill/NEON_ASVs/ITS/'

#Function will loop over multiple sites and dates. These should be taken from a nested site*date list.
site <- 'SCBI'
date <- '2014-09'

#these IDs need to come from the data product. This vector is based on a unqiue site*date combination.
mgm.id <- c('mgm4788103.3','mgm4788105.3')

#make site-date directory if it doesnt already exist.
#if upstream site directory doesnt exist, it will also be created.
cur.dir <- paste0(main.path,site,'/',date)
pre <- paste0('mkdir -p ', cur.dir)
system(pre)

#output path.
output.path <- paste0(cur.dir,'/ASV_',site,'.',date,'.rds')
#check if output file already present. If it is, stop doing this, and report this.
#turn this on once you running a function loop.
#if(file.exists(output.path) == T){
#  cat(paste0(site,' ',date,' ASV table already present. you good! moving on.\n'))
#  next
#  }

#download sequences from MG-RAST. These reads are paired, and presumably quality filtered.
cat(paste0('Downloading sequences for ',site,' ',date,'...\n'))
for(i in 1:length(mgm.id)){
  mgrast.link <- paste0('http://api.metagenomics.anl.gov/1/download/',mgm.id[i],'?file=050.1')
  file <- paste0(cur.dir,'/',mgm.id[i],'.fasta')
  pre <- paste0('curl ',mgrast.link,' > ',file)
  system(pre)
}
cat('Download complete.')

#get names of all sequence files in cur.dir.
files <- sort(list.files(cur.dir, pattern = '.fasta', full.names = T))
sample.names <- sapply(strsplit(basename(files), ".fasta"), `[`, 1)


#Build an ASV table from all sequence files.
cat(paste0('Building ASV table for ',site,' ',date,'...'))
for(i in 1:length(sample.names)){
  #get just the sequences into a separate file.
    file <- paste0(cur.dir,'/'    ,sample.names[i],'.fasta')
  s.file <- paste0(cur.dir,'/seq.',sample.names[i],'.fasta')
  pre <- paste0("sed -n '0~2p' ",file,' > ',s.file)
  system(pre)
  
  #Convert sequences to a sequence table. Write table to a file. No need to load all sequences into R memory.
  c.file <- paste0(cur.dir,'/counts.',sample.names[i],'.txt')   #sample-specific ASV table.
  pre <- paste0('cat ',s.file,' | sort | uniq -c > ', c.file)
  system(pre)
  asv <- data.table::data.table(read.table(c.file))
  colnames(asv) <- c(sample.names[i],'seq')
  data.table::setkey(asv,seq)

  #merge multiple sample-specific ASV tables.
  if(i == 1){ out = asv}
  if(i > 1){ out <- merge(out, asv, all = T)} #this merge really requires data.table be loaded into the environment.
  
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
colnames(t.out) <- out$seq

#convert from numeric dataframe to integer matrix. this is important for dada2 commands downstream.
t.out <- as.matrix(t.out)
t.out <- apply (t.out, c (1, 2), function (x) {(as.integer(x))})

#use dada2 to remove chimeras. We are going to do this in a different function.
#need to do the same for assigning RDP taxonomy in dada2.
#It takes forever and its doing redundant stuff if you merge multiple tables, which you will.
#This takes a while. Could be sped up with more cores I think.
#cat('Removing chimeric sequences (this takes a while)...')
#t.out <- dada2::removeBimeraDenovo(t.out, method="consensus", multithread=TRUE, verbose=TRUE)
#cat('Chimeric sequences removed!')

#save ASV table, remove downloaded sequence data.
saveRDS(t.out,output.path)
system(paste0('rm ',cur.dir,'/*.fasta'))

cat(paste0(site,' ',date,' ASV table constructed.'))


#This is where function will end. Wrap this in date and site loops.
#After this function we need a second that assigns taxonomy and writes a file.
#A third that assigns functional groups and writes a file.
#A fourth that generates OTU tables based on taxonomy and functional types.