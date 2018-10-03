#trimming primers from fastq files for bahram 2018.
#relies in bbduk from the bbmap package.
#This only works on scc, as java is not installed or loaded or up to date or something on pecan2.
#Note- you only retain ALL reads if you sent min_len=1. Default is min_len=10, which drops a bunch (though, we probably do want to discard these).
rm(list=ls())
source('paths.r')

#Input Primers: these include degenerate (variable) bases, according to IUPAC standard. bbduk can handle this.
fwd.primer <- 'GTGYCAGCMGCCGCGGTAA'
rev.primer <- 'GGACTACNVGGGTWTCTAAT'

#Input path to sequences:
seq.dir <- bahram.seq.dir


#once directory and primers specified, below code is independent on user supplied parameters.
files <- list.files(seq.dir)
files <- files[grep('.fastq',files)]
#subset forward and reverse reads.
fwd.files <- files[grep('_1.fastq',files)]
rev.files <- files[grep('_2.fastq',files)]

#path to bbduk tool w/in the bbmap directory.
bbduk.path <- 'NEFI_functions/bbmap/bbduk.sh'

#trim left, "forward trim".
for(i in 1:length(fwd.files)){
  sample.name1 <- fwd.files[i]
  sample.name2 <- rev.files[i]
  input_sample.path1 <- paste0(seq.dir,sample.name1)
  input_sample.path2 <- paste0(seq.dir,sample.name2)
  output_sample.path1 <- paste0(seq.dir,'q.trim.L/',sample.name1)
  output_sample.path2 <- paste0(seq.dir,'q.trim.L/',sample.name2)
  cmd <- paste0(bbduk.path, 
                ' in1=',input_sample.path1,
                ' in2=',input_sample.path2,
                ' out1=',output_sample.path1,
                ' out2=',output_sample.path2,
                ' literal=',fwd.primer,
                ' ktrim=l k=10 ordered=t minlen=1')
  system(cmd)  
}
#trim right, "reverse trim". bbduk automatically looks for reverse complement.
for(i in 1:length(fwd.files)){
  sample.name1 <- fwd.files[i]
  sample.name2 <- rev.files[i]
  input_sample.path1 <- paste0(seq.dir,'q.trim.L/',sample.name1)
  input_sample.path2 <- paste0(seq.dir,'q.trim.L/',sample.name2)
  output_sample.path1 <- paste0(seq.dir,'q.trim.R/',sample.name1)
  output_sample.path2 <- paste0(seq.dir,'q.trim.R/',sample.name2)
  cmd <- paste0(bbduk.path,
                ' in1=',input_sample.path1,
                ' in2=',input_sample.path2,
                ' out1=',output_sample.path1,
                ' out2=',output_sample.path2,
                ' literal=',rev.primer,
                ' ktrim=l k=10 ordered=t minlen=1')
  system(cmd)  
}

#clean up.
cmd <- paste0('rm -rf ',seq.dir,'q.trim.L')
system(cmd)
cmd <- paste0('mv ',seq.dir,'q.trim.R ',seq.dir,'q.trim')
system(cmd)

