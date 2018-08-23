#assign taxonomy to tedersoo sequences.

#clear environment, load packages, functions and paths.----
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')

#load ASV table, set output path.----
#this needs a lot of memory.
d <- readRDS(ted_2014_SV.table.path)
output.path <- ted_2014_tax.path

#download unite training set.----
cat('downloading UNITE database...\n')
unite_url <- 'https://files.plutof.ut.ee/doi/B2/07/B2079372C79891519EF815160D4467BBF4AF1288A23E135E666BABF2C5779767.zip'
unite_path.zip <- paste0(data.dir,'unite.fa.zip')
unite_path     <- paste0(data.dir,'sh_general_release_dynamic_01.12.2017.fasta')
cmd <- paste0('curl ',unite_url,' > ',unite_path.zip)
system(cmd)
cmd <- paste0('unzip ',unite_path.zip)
system(cmd)
cat('UNITE download complete.\n')

#assign taxonomy.----
tax <- dada2::assignTaxonomy(colnames(d),unite_path, multithread = T)

#save output.----
saveRDS(tax, output.path)

#cleanup.----
system('rm ',unite_path.zip)
system('rm ',unite_path)

#end script.
