# assign functional groups to 16S taxa from NEON.

rm(list=ls())
library(data.table)
source('paths.r')

# read in csv with classifications. 
fg <- read.csv(paste0(pecan_gen_16S_dir, "bacteria_func_groups.csv"))

# load Bahram SV table as otu file
otu <- readRDS(NEON_dada2_SV_table.path)

# load Bahram taxonomy
tax <- readRDS(NEON_dada2_tax_table.path)

# remove leading "k__" in taxonomy.
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

# for column names to be lower case.
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'Bacteria'|tax$kingdom == 'Archaea',] # only removes ~500 counts
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]


# start with copiotroph/oligotroph.
c_o_groups <- fg[fg$Classification.system=="Copiotroph_oligotroph",]
copiotrophs <- c_o_groups[c_o_groups$Classification=="Copiotroph",]$Taxon
oligotrophs <- c_o_groups[c_o_groups$Classification=="Oligotroph",]$Taxon
tax$group <- NA

# first assign at phylum level
tax[which(tax$phylum %in% copiotrophs),]$group <- "copiotroph"
tax[which(tax$phylum %in% oligotrophs),]$group <- "oligotroph"

# then class level
tax[which(tax$class %in% copiotrophs),]$group <- "copiotroph"
tax[which(tax$class %in% oligotrophs),]$group <- "oligotroph"

# then order level
tax[which(tax$order %in% copiotrophs),]$group <- "copiotroph"
tax[which(tax$order %in% oligotrophs),]$group <- "oligotroph"

# then family level
tax[which(tax$family %in% copiotrophs),]$group <- "copiotroph"
tax[which(tax$family %in% oligotrophs),]$group <- "oligotroph"

# then genus level
tax[which(tax$genus %in% copiotrophs),]$group <- "copiotroph"
tax[which(tax$genus %in% oligotrophs),]$group <- "oligotroph"


#Get seq abundances of copiotrophs vs oligotrophs.----
classification <- c("copiotroph", "oligotroph")
cop_olig <- list()
k <- data.table(cbind(tax,t(otu)))
for(i in 1:length(classification)){
  z <- k[group == classification[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  cop_olig[[i]] <- out
}
cop_olig <- data.frame(t(do.call('rbind',cop_olig)))
colnames(cop_olig) <- classification
seq_total <- colSums(k[,start:ncol(k)])
other <- seq_total - rowSums(cop_olig)
cop_olig <- cbind(other,cop_olig)
cop_olig <- list(cop_olig,seq_total)
names(cop_olig) <- c('abundances','seq_total')
cop_olig$rel.abundances <- cop_olig$abundances / cop_olig$seq_total
saveRDS(cop_olig, NEON_cop_olig_16S.path)
