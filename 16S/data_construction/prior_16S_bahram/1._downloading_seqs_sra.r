#ripping Bahram et al. 2018 Nature 16S sequences from the SRA database.
#note: this could go faster using the aspera connect tools. But this works.
#clear environment, load packages.
rm(list=ls())
source('paths.r')
#source('http://bioconductor.org/biocLite.R') #you need to install the package this way if you don't already have it.
#biocLite('SRAdb')
library(SRAdb)
srafile = getSRAdbFile()
con = dbConnect(RSQLite::SQLite(),srafile) #this takes a frickin minute. 2.3G.

#Grab ftp addresses and acession numbers in case you want them later.
d <- listSRAfile('ERP021922', con)
#write.csv(d,'Bahram_2018_16S_accessions.csv')

#Download all data to bahram sequence directory. Path is specified in paths.r
#this is like 10G.
getSRAfile('ERP021922', con, fileType = 'fastq',destDir = bahram.seq.dir)

#unzip all fastq.gz files.
cmd <- paste0('gunzip ',bahram.seq.dir,'/*.gz')
system(cmd)

#clean up. Remove SRAdb sqlite file.
system('rm -f SRAmetadb.sqlite')
system('rm -f SRAmetadb.sqlite.gz')

#end script.
"$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"