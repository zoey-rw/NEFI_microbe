# assign functional groups to 16S taxa from Bahram.
# 1. prep data.
# 2. assign copiotroph/oligotroph groups from lit review (csv file, downloaded from google sheet)
# 3. assign N-cycling groups from lit review, and then from Albright et al. 2018 dataset.
# 4. assign C-cycling groups from lit review.

rm(list=ls())
library(data.table)
library(readxl)
library(tidyr)
library(stringr)
source('paths.r')

########## 1. prep data. ############

# load csv with literature-review classifications. 
fg <- read.csv(paste0(pecan_gen_16S_dir, "bacteria_func_groups.csv"))

# load csv from Albright with N-cycle pathway presence/absence.
N_cyclers <- read_excel(paste0(pecan_gen_16S_dir, "Npathways_Albright2018.xlsx"))

# load csv from Berlemont and Martiny with cellulolytic pathway presence/absence
cell <- read.csv(paste0(pecan_gen_16S_dir, "cellulolytic_Berlemont.csv"))

# load Bahram SV table as otu file
otu <- readRDS(bahram_dada2_SV_table.path)

# load Bahram taxonomy
tax <- readRDS(bahram_dada2_tax_table.path)

# load metadata
metadata <- readRDS(bahram_metadata.path)

# remove leading "k__" in taxonomy.
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

# subset otu table and tax table to only include observations in map file
metadata$Run <- as.character(metadata$Run)
otu <- otu[rownames(otu) %in% metadata$Run,]
metadata <- metadata[metadata$Run %in% rownames(otu),]
# order OTU table to match the mapping file
otu <- otu[order(rownames(otu), metadata$Run),]

# for column names to be lower case.
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'Bacteria'|tax$kingdom == 'Archaea',] # only removes ~500 counts
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]
tax_save <- tax # just so we have this taxonomic table for later.







########## 2. assign copiotroph/oligotroph groups. ############

# read in groups from csv
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

#Get seq abundances of copiotrophs vs oligotrophs, in one dataframe.----
classification <- c("copiotroph", "oligotroph")
cop_olig <- list()
k <- data.table(cbind(tax,t(otu)))
for(i in 1:length(classification)){
  z <- k[group == classification[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  cop_olig[[i]] <- out
}
cop_olig_save <- cop_olig
cop_olig <- data.frame(t(do.call('rbind',cop_olig)))
colnames(cop_olig) <- classification
seq_total <- colSums(k[,start:ncol(k)])
other <- seq_total - rowSums(cop_olig)
cop_olig <- cbind(other,cop_olig)
cop_olig <- list(cop_olig,seq_total)
names(cop_olig) <- c('abundances','seq_total')
cop_olig$rel.abundances <- cop_olig$abundances / cop_olig$seq_total
saveRDS(cop_olig, prior_cop_olig_abundances.path)


########## 3. assign nitrogen-cycling groups. ############

# Set up N-cycle dataset from Albright 2008

tax <- tax_save 

# Remove all columns except for taxonomy, environment, and pathways
N_cyclers <- N_cyclers[,-c(1:2,4,11:14,16:17)] 

# Rename some pathways
setnames(N_cyclers, 
         c("Nitrogen Fixation", "Assimilatory Nitrite to ammonia", "Dissimilatory Nitrite to Ammonia", "Assimilatory Nitrate to Nitrite","Dissimilatory Nitrate to Nitrite"), 
         c("N_fixation", "Assim_nitrite_reduction", "Dissim_nitrite_reduction", "Assim_nitrate_reduction","Dissim_nitrate_reduction"))

# Treat "incomplete" pathways as if they are absent. 
N_cyclers[N_cyclers=="complete"] <- 1
N_cyclers[N_cyclers=="incomplete" | N_cyclers=="None"] <- 0

# Grouping partial nitrification pathway with nitrification, and partial denitrification with denitrification.
N_cyclers[N_cyclers$Partial_Nitrification==1,]$Nitrification <- 1
N_cyclers[N_cyclers$Partial_NO==1 | N_cyclers$Partial_N2O ==1 | N_cyclers$Partial_N2==1,]$Denitrification <- 1

# Now we can remove the "partial" columns.
N_cyclers[,c("Partial_Nitrification", "Partial_NO", "Partial_N2O", "Partial_N2")] <- NULL

# # Get all possible combinations of the core seven pathways - just curious.
# possible_combos <- expand(N_cyclers, Assim_nitrite_reduction, Dissim_nitrite_reduction, Assim_nitrate_reduction, Dissim_nitrate_reduction, N_fixation, Nitrification, Denitrification)
# # nrow(possible_combos)
# # [1] 128
# actual_combos <- expand(N_cyclers, nesting(Assim_nitrite_reduction, Dissim_nitrite_reduction, Assim_nitrate_reduction, Dissim_nitrate_reduction, N_fixation, Nitrification, Denitrification))
# # nrow(actual_combos)
# # [1] 50


# For Albright dataset:
# get value by genus; if pathway is present in more than half of species, it is present for that genus
pathway_names <- colnames(N_cyclers[9:15])
# N_cycle_genera <- data.frame(matrix(ncol=0,nrow=0))
# for (g in 1:length(unique(N_cyclers$Genus))) {
#   i <- unique(N_cyclers$Genus)[g]
#   species <- N_cyclers[N_cyclers$Genus==i,]
#   nspecies <- nrow(species)
#   out <- data.frame(Genus = i,
#                     nrow(species[species$Assim_nitrite_reduction == 1,])/nspecies,
#                     nrow(species[species$Dissim_nitrite_reduction == 1,])/nspecies,
#                     nrow(species[species$Assim_nitrate_reduction == 1,])/nspecies,
#                     nrow(species[species$N_fixation == 1,])/nspecies,
#                     nrow(species[species$Dissim_nitrate_reduction == 1,])/nspecies,
#                     nrow(species[species$Nitrification == 1,])/nspecies,
#                     nrow(species[species$Denitrification == 1,])/nspecies
#   )
#   colnames(out) <- c("Genus", pathway_names)
#   out$Genus <- as.character(i)
#   out[out >= .5] <- 1
#   out[out < .5] <- 0
#   out$Genus <- as.character(i)
#   N_cycle_genera <- rbind(N_cycle_genera, out)
# }

##### Assign taxa to functional groups ##### 

# create pathway columns
tax[,pathway_names] <- NA

# check if sample genus is in classification data, and that a classified pathway is present; 
# assign those genera a present pathway
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  # Classifications from literature search (multiple taxon levels)
  has_pathway <- fg[fg$Classification==p,]$Taxon
  
  if (nrow(tax[tax$phylum %in% has_pathway,]) != 0){
    tax[tax$phylum %in% has_pathway,][p] <- 1
  }
  if (nrow(tax[tax$class %in% has_pathway,]) != 0){
    tax[tax$class %in% has_pathway,][p] <- 1
  }
  if (nrow(tax[tax$order %in% has_pathway,]) != 0){
    tax[tax$order %in% has_pathway,][p] <- 1
  }
  if (nrow(tax[tax$family %in% has_pathway,]) != 0){
    tax[tax$family %in% has_pathway,][p] <- 1
  }
  if (nrow(tax[tax$genus %in% has_pathway,]) != 0){
    tax[tax$genus %in% has_pathway,][p] <- 1
  }
  # genus + species must match any full species name
  if (nrow(tax[which(paste(tax$genus, tax$species) %in% has_pathway),]) != 0){
    tax[which(paste(tax$genus, tax$species) %in% has_pathway),][p] <- 1
  }
  # Classifications from Albright et al. 2018 dataset (Genus-level only; reduced from species-level)
  has_pathway <- N_cyclers[N_cyclers[p] == 1,]$Genus # not using .5 cutoff - one species with pathway is enough to classify genus

 if (nrow(tax[which(tax$genus %in% has_pathway),][p]) != 0){
    tax[which(tax$genus %in% has_pathway),][p] <- 1
  }
}

# check how many ended up with classifications.
tax_classified <- tax[,8:14]
no_pathways <- tax_classified[apply(tax_classified, 1, function(x) !any(x == 1)),]
nrow(no_pathways)
# 133599 # when all genera are included, not .5 cutoff
some_pathway <- tax_classified[apply(tax_classified, 1, function(x) any(x == 1)),]
nrow(some_pathway)
# 22643 # not .5 cutoff


#Get seq abundances of each pathway
allpathways <- list()
k <- data.table(cbind(tax,t(otu)))
for(i in 1:length(pathway_names)){
  pathways <- list()
  n <- pathway_names[i]
  z <- k[k[[n]] == 1,]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  pathways[[i]] <- out
  pathways <- t(do.call('rbind',pathways))
  colnames(pathways) <- pathway_names[i]
  seq_total <- colSums(k[,start:ncol(k)])
  other <- seq_total - rowSums(pathways)
  pathways <- cbind(other,pathways)
  pathways <- list(pathways,seq_total)
  names(pathways) <- c('abundances','seq_total')
  pathways$rel.abundances <- pathways$abundances / pathways$seq_total
  allpathways[[i]] <- pathways
}
#saveRDS(allpathways, prior_N_cyclers_abundances.path)
saveRDS(allpathways, prior_N_cyclers_abundances.path)





########## 4. assign C-cycling groups. ###########

# format dataset of cellulolytic taxa from Berlemont et al.
rownames(cell) <- cell$Strain

# these have enzymes for beta-glucosidase, to use products from cellulysis - opportunists
bg <- cell[,c(3,4)]
bg <- bg[apply(bg, 1, function(x) any(x == 1)),]

# taxa with any enzymes for cellulolysis 
cellulolytic <- cell[,c(5:7,9:13)]
cellulolytic$is.cell <- 0
cellulolytic[apply(cellulolytic, 1, function(x) any(x == 1)),]$is.cell <- 1
cellulolytic$genus <- word(rownames(cellulolytic), 1)
cellulolytic[cellulolytic$genus == "Candidatus",]$genus <- 
  word(rownames(cellulolytic[cellulolytic$genus == "Candidatus",]), 1, 2)

cellulolytic_genera <- data.frame(matrix(ncol=0,nrow=0))
for (g in 1:length(unique(cellulolytic$genus))) {
  i <- unique(cellulolytic$genus)[g]
  species <- cellulolytic[cellulolytic$genus==i,]
  nspecies <- nrow(species)
  out <- data.frame(Genus = i,
                    nrow(species[species$is.cell == 1,])/nspecies
  )
  colnames(out) <- c("Genus", "is.cell")
  out$Genus <- as.character(i)
  out[out >= .5] <- 1
  out[out < .5] <- 0
  out$Genus <- as.character(i)
  cellulolytic_genera <- rbind(cellulolytic_genera, out)
}


########## assign carbon-cycling taxa ############

tax <- tax_save
#pathway_names <- unique(fg[fg$Classification.system=="C cycling",]$Classification) # don't want all right now
pathway_names <- c("Cellulolytic", "Chitinolytic", "Lignolytic", "Methanotroph") # removed "Methanogen" - only present in one sample

tax[,pathway_names] <- NA

# check if sample genus is in classification data, and that a classified pathway is present; 
# assign those genera a present pathway
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  print(p)
  
  # Classifications from literature search (multiple taxon levels)
  has_pathway <- fg[which(fg$Classification==p),]$Taxon
  
  if (nrow(tax[tax$phylum %in% has_pathway,]) != 0){
    tax[tax$phylum %in% has_pathway,][p] <- 1
  }
  if (nrow(tax[tax$class %in% has_pathway,]) != 0){
    tax[tax$class %in% has_pathway,][p] <- 1
  }
  if (nrow(tax[tax$order %in% has_pathway,]) != 0){
    tax[tax$order %in% has_pathway,][p] <- 1
  }
  if (nrow(tax[tax$family %in% has_pathway,]) != 0){
    tax[tax$family %in% has_pathway,][p] <- 1
  }
  if (nrow(tax[tax$genus %in% has_pathway,]) != 0){
    tax[tax$genus %in% has_pathway,][p] <- 1
  } 
  # genus + species must match any full species name
  if (nrow(tax[which(paste(tax$genus, tax$species) %in% has_pathway),]) != 0){
    tax[which(paste(tax$genus, tax$species) %in% has_pathway),][p] <- 1
  }
  # Classifications from Berlemont et al. 2018 dataset (Genus-level only; reduced from species-level)
  if (p == "Cellulolytic") {
    has_pathway <- cellulolytic[cellulolytic$is.cell == 1,]$genus # not using .5 cutoff - one species with pathway is enough to classify genus
  if (nrow(tax[which(tax$genus %in% has_pathway),][p]) != 0){
    tax[which(tax$genus %in% has_pathway),][p] <- 1
  }
  } # close cellulolytic section
  
}

# # check how many ended up with classifications.
tax_classified <- tax[,8:11]
no_pathways <- tax_classified[apply(tax_classified, 1, function(x) !any(x == 1)),]
nrow(no_pathways)
# [1] 151048 #no cutoff
some_pathway <- tax_classified[apply(tax_classified, 1, function(x) any(x == 1)),]
nrow(some_pathway)
# [1] 5194 #no cutoff


#Get seq abundances of each pathway
allpathways <- list()
k <- data.table(cbind(tax,t(otu)))
for(i in 1:length(pathway_names)){
  pathways <- list()
  n <- pathway_names[i]
  z <- k[k[[n]] == 1,]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  pathways[[i]] <- out
  pathways <- t(do.call('rbind',pathways))
  colnames(pathways) <- pathway_names[i]
  seq_total <- colSums(k[,start:ncol(k)])
  other <- seq_total - rowSums(pathways)
  pathways <- cbind(other,pathways)
  pathways <- list(pathways,seq_total)
  names(pathways) <- c('abundances','seq_total')
  pathways$rel.abundances <- pathways$abundances / pathways$seq_total
  allpathways[[i]] <- pathways
}
saveRDS(allpathways, prior_C_cyclers_abundances.path)


