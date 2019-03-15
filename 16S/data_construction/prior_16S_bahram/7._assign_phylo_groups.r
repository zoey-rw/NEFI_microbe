#Aggregate cosmpolitan taxa in Bahram.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
#source('NEFI_functions/common_group_quantification.r')

library(RCurl)
# source function from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/common_group_quantification.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

library(data.table)

#set output paths.----
output.path <- bahram_16S_common_phylo_groups_list.path

#load data.----
map <- readRDS(bahram_metadata.path)
otu <- readRDS(bahram_dada2_SV_table_rare.path)
tax <- readRDS(bahram_dada2_tax_table.path)

#subset otu file to match mapping file.----
otu <- otu[rownames(otu) %in% map$Run,]

#format taxonomy table.----
#remove leading "k__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

#for column names to be lower case.
colnames(tax) <- tolower(colnames(tax))

#remove taxa that do not assign to fungi from tax and otu table.----
tax <- as.data.frame(tax)
tax <- tax[tax$kingdom == 'Bacteria'|tax$kingdom == 'Archaea',] # only removes ~500 counts
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

#get each level of taxonomy output.----
phylum.out <- common_group_quantification(otu,tax,unique(tax$phylum),'phylum', samp_freq = .95)
class.out <- common_group_quantification(otu,tax,unique(tax$class),'class', samp_freq = .95 )
order.out <- common_group_quantification(otu,tax,unique(tax$order),'order', samp_freq = .95)
family.out <- common_group_quantification(otu,tax,unique(tax$family),'family', samp_freq = .95)
genus.out <- common_group_quantification(otu,tax,unique(tax$genus),'genus', samp_freq = .95)

#save output.----
all_taxa_out <- list(phylum.out,class.out,order.out,family.out,genus.out)
names(all_taxa_out) <- c('phylum','class','order','family','genus')
saveRDS(all_taxa_out,output.path) 
