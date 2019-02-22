#NEON taxonomic group tables for all phylogenetic levels - 16S.
#common = found in greater than 50% of prior samples.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')
#source('NEFI_functions/common_group_quantification.r')

# source hierarch means function
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/common_group_quantification.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#set output path.----
output.path <- NEON_16S_phylo_groups_abundances.path

#load data.----
sv <- readRDS(NEON_dada2_SV_table.path)
map <- readRDS(hierarch_filled_16S.path)
map <- map$core.obs
tax <- readRDS(NEON_dada2_tax_table.path)

#get reference common phyla from prior.----
ref <- readRDS(bahram_16S_common_phylo_groups_list.path)

#check that things are properly ordered.----
if(sum(rownames(tax) == colnames(sv)) != ncol(sv)){
  cat('Stop. rownames of tax-fun table do not match column names of SV table.\n')
}

# prep taxonomic table
# remove leading "k__" in taxonomy.
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

# for column names to be lower case.
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))


#aggregate at multiple phylogenetic levels.----
levels <- c('phylum','class','order','family','genus')
output <- list()
for(i in 1:length(levels)){
  level <- levels[i]
  reference <- data.table(ref[[i]]$group_frequencies)
  reference <- as.character(reference[sample_frequency > 0.5]$groups)
  output[[i]] <- common_group_quantification(sv, tax, reference, level, ref_filter = T)
}
names(output) <- names(ref)

#save output.----
saveRDS(output, output.path)
