#remove chimeras with dada2.
#1. clear environment, load functions and packages.----
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(data.table)

#2. Setup paths.----
#get sequence file paths
seq.path <- NEON_ITS.dir
SV_pre.chimera.path <- paste0(seq.path,'SV_pre.chimera_table.rds')
#output file paths.
output_filepath1 <-  paste0(seq.path,'SV_table.rds')
output_filepath2 <- NEON_SV.table.path

#3. load data, remove chimeras.----
t.out <- readRDS(SV_pre.chimera.path)
cat('Removing chimeras...\n')
tic()
t.out_nochim <- dada2::removeBimeraDenovo(t.out, method = 'consensus', multithread = TRUE)
cat('Chimeras removed.\n')
toc()

#4. Final save and cleanup.----
#save output.
cat('Saving output.../n')
tic()
saveRDS(t.out_nochim, output_filepath1)
saveRDS(t.out_nochim, output_filepath2)
toc()
cat('Output saved./n')

#clean up if t.out_nochim was created.
if(exists('t.out_nochim')==T){
  cmd <- paste0('rm -f ',SV_pre.chimera.path)
  system(cmd)
}

#end script.
