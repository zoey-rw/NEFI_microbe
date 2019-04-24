# assign phylogenetic lineages for major groups

rm(list=ls())
source('paths.r')
source('NEFI_functions/common_group_quantification.r')
library(data.table)

#set output paths.----
output.path <- bahram_16S_phylo_lineage_abun_16S.path

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

dim(tax[tax$phylum %in% c("Acidobacteria"),])
acido <- tax[tax$phylum %in% c("Acidobacteria"),]
rownames(acido) <- NULL
acido <- acido[,!colnames(acido) %in% c("species")]
acido <- acido[!duplicated(acido),]
acido <- acido[!is.na(acido$class) & acido$class != "",]
acido[is.na(acido)] <- ""
acido <- acido[!duplicated(acido),]

library(data.tree)
#acido$R2 <- 
# site level R2: 
  # acido$phylum "Acidobacteria"  .4, 0% in predictive interval
  # acido$class "Acidobacteria-6"  .71, 91% in predictive interval
  # acido$class "iii1-8" .34, 100% in predictive interval
  # acido$order "iii1-15" .66, 91% in predictive interval
  # acido$order "RB41" .73 82%
  # acido$class "Solibacteres" .79 100%
  # acido$order "Solibacterales" .75 100%
  # 
acido$pathString <- paste(acido$phylum, acido$class, acido$order, acido$family, acido$genus, sep="/")
acido.tree <- as.Node(acido)
print(acido.tree)
plot(as.dendrogram(acido.tree, edgetext = T))

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