#Downloading raw ITS sequences from NEON.
#This needs to be run on the SCC via the qsub script so it will properly load python and qiime packages.
#this code has been tested line by line, but the qsub script on the scc has not.
#clear environment, source paths and packages.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')

#load sheet with links to data from L. Stanish.
d <- read.csv(NEON_ITS_link_file)

#subset to only ITS reads.
d <- d[grep('ITS',d$rawDataFilePath),]

#grab unique URLs. Only two.
neon_paths <- as.character(unique(d$rawDataFilePath))

#download.
cat('Downloading sequences...\n')
tic()
for(i in 1:length(neon_paths)){
  path <- neon_paths[i]
  name <- sub('.*.neonscience.org/','',path)
  name <- sub('.*/','',name)
  output.path <- paste0(NEON_ITS.dir,name)
  cmd <- paste0('wget -O ',output.path,' ',path, '  --no-check-certificate')
  #cmd <- paste0('wget ',path,' -O ',NEON_ITS.dir)
  #cmd <- paste0('wget wget -e robots=off -r --no-parent ',NEON_ITS.dir,' ',path,' --no-check-certificate')
  system(cmd)
  #unzip.
  cmd <- paste0('unzip ',output.path,' -d ',NEON_ITS.dir)
  system(cmd)
  #cleanup zip file.
  cmd <- paste0('rm ',output.path)
  system(cmd)
}
cat('Download complete! ')
toc()

#replace all spaces in file names.
cmd <- paste0('rename " " "_" ',NEON_ITS.dir,'*')
system(cmd)

#break into per-sample sequence files so I can hack this in parallel.
files <- list.files(NEON_ITS.dir)
files <- files[grep('.fna',files)]
for(i in 1:length(files)){
  this.file <- paste0(NEON_ITS.dir,files[i])
  cmd <- paste0('split_sequence_file_on_sample_ids.py -i ',this.file,' -o ',NEON_ITS.dir)
  system(cmd)
}

#clean up raw files
for(i in 1:length(files)){
  this.file <- paste0(NEON_ITS.dir,files[i])
  cmd <- paste0('rm ',this.file)
  system(cmd)
}

#end script.
