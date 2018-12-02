#Aggregate cosmpolitan genera.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/common_group_quantification.r')
library(data.table)

#set output paths.----
#cosmo_output.path <- tedersoo_ITS_cosmo_genera_list.path
#phyla_output.path <- tedersoo_ITS_phyla_list.path
output.path <- tedersoo_ITS_common_phylo_groups_list.path

#load data.----
map <- readRDS(tedersoo_ITS_clean_map.path)
otu <- readRDS(ted_2014_SV.table.path)
tax <- readRDS(ted_2014_tax.path)

#subset otu file to match mapping file.----
otu <- otu[rownames(otu) %in% map$SRR.id,]

#format taxonomy table.----
#remove leading "k__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

#for column names to be lower case.
colnames(tax) <- tolower(colnames(tax))

#remove taxa that do not assign to fungi from tax and otu table.----
tax <- tax[tax$kingdom == 'Fungi',]
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

#get each level of taxonomy output.----
phylum.out <- common_group_quantification(otu,tax,unique(tax$phylum),'phylum')
 class.out <- common_group_quantification(otu,tax,unique(tax$class ),'class' )
 order.out <- common_group_quantification(otu,tax,unique(tax$order ),'order' )
family.out <- common_group_quantification(otu,tax,unique(tax$family),'family')
 genus.out <- common_group_quantification(otu,tax,unique(tax$genus ),'genus' )

#save output.----
all_taxa_out <- list(phylum.out,class.out,order.out,family.out,genus.out)
names(all_taxa_out) <- c('phylum','class','order','family','genus')
saveRDS(all_taxa_out,output.path) 

#saveRDS( genus.out, cosmo_output.path)
#saveRDS(phylum.out, phyla_output.path)
