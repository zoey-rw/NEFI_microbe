#NEON taxonomic group tables for all phylogenetic levels - 16S.
#common = found in greater than 50% of prior samples.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')
source('NEFI_functions/common_group_quantification.r')

# source function
# library(RCurl)
# script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/common_group_quantification.r", ssl.verifypeer = FALSE)
# eval(parse(text = script))

#set output path.----
output.path <- NEON_16S_phylo_fg_abundances.path

#load data.----
otu <- readRDS(NEON_dada2_SV_table_rare.path)
tax <- readRDS(NEON_tax_fg_16S.path)

#get reference common phyla from prior.----
ref <- readRDS(bahram_16S_common_phylo_fg_abun.path)

#check that things are properly ordered.----
if(sum(rownames(tax) == colnames(otu)) != ncol(otu)){
  cat('Stop. rownames of tax-fun table do not match column names of SV table.\n')
}

#get each level of taxonomy output.----
of_interest <- c('phylum','class','order','family','genus',colnames(tax)[8:19])
all_taxa_out <- list()
for(i in 1:length(of_interest)){
  reference <- data.table(ref[[i]]$group_frequencies)
  reference <- as.character(reference[sample_frequency > 0.95]$groups)
  all_taxa_out[[i]] <- common_group_quantification(otu,
                                                   tax,
                                                   reference,
                                                   of_interest[i],
                                                   ref_filter = T)
}
names(all_taxa_out) <- of_interest

#save output.----
saveRDS(all_taxa_out,output.path) 
