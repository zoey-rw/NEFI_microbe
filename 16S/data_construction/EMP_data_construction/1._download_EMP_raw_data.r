#Download EMP deblurred ESV table and mapping file.
#clear environment, load packages and source file paths.
rm(list=ls())
source('paths.r')

#download EMP ESV table.
emp_esv.url <- 'ftp://ftp.microbio.me/emp/release1/otu_tables/deblur/emp_deblur_150bp.release1.biom'
command <- paste0('curl -o ',emp_esv.path,' ',emp_esv.url)
system(command)

#download EMP mapping file.
emp_map.url <- 'ftp://ftp.microbio.me/emp/release1/mapping_files/emp_qiime_mapping_release1.tsv'
command <- paste0('curl -o ',emp_map.path,' ',emp_map.url)
system(command)
