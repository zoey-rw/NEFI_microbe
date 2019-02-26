# assign functional groups to 16S taxa from NEON.
# 1. prep data.
# 2. assign copiotroph/oligotroph groups from lit review (csv file, downloaded from google sheet)
# 3. assign N-cycling groups from lit review, and then from Albright et al. 2018 dataset.
# 4. assign C-cycling groups from lit review.

# clear environment, load libraries.
rm(list=ls())
library(data.table)
library(readxl)
library(tidyr)
library(stringr)
source('paths.r')

# read in reference and taxonomic data.

# load csv with classifications
fg <- read.csv(paste0(pecan_gen_16S_dir, "bacteria_func_groups.csv"))
# load excel file from Albright with N-cycle pathway presence/absence
N_cyclers <- read_excel(paste0(pecan_gen_16S_dir, "Npathways_Albright2018.xlsx"))
# load csv from Berlemont and Martiny with cellulolytic pathway presence/absence
cell <- read.csv(paste0(pecan_gen_16S_dir, "cellulolytic_Berlemont.csv"))
# load NEON SV table as otu file
otu <- readRDS(NEON_dada2_SV_table_rare.path)
# load NEON taxonomy
tax <- readRDS(NEON_dada2_tax_table.path)


########## 1. prep taxonomic data. ############

# remove leading "k__" in taxonomy.
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

# for column names to be lower case.
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'Bacteria'|tax$kingdom == 'Archaea',] # only removes ~25 counts
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]
tax_save <- tax # just so we have this taxonomic table for later.



########## 2. assign copiotroph/oligotroph groups. ############
tax <- tax_save
tax$group <- NA

# read in groups from csv
c_o_groups  <-  fg[fg$Classification.system == "Copiotroph_oligotroph",]
copiotrophs <-  c_o_groups[c_o_groups$Classification == "Copiotroph",]$Taxon
oligotrophs <-  c_o_groups[c_o_groups$Classification == "Oligotroph",]$Taxon

# taxon assignments
levels <- c("phylum", "class", "order", "family", "genus")
for (i in 1:length(levels)) {
  p <- levels[i]
  tax[which(tax[[p]] %in% copiotrophs),]$group <- "copiotroph"
  tax[which(tax[[p]] %in% oligotrophs),]$group <- "oligotroph"
}

#Get seq abundances of copiotrophs vs oligotrophs, in one dataframe.----
classification <- c("copiotroph", "oligotroph")
cop_olig <- list()
k <- data.table(cbind(tax, t(otu)))
for (i in 1:length(classification)) {
  z <- k[group == classification[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[, start:ncol(z)])
  cop_olig[[i]] <- out
}
cop_olig <- data.frame(t(do.call('rbind', cop_olig)))
colnames(cop_olig) <- classification
seq_total <- colSums(k[, start:ncol(k)])
other <- seq_total - rowSums(cop_olig)
cop_olig <- cbind(other, cop_olig)
cop_olig <- list(cop_olig, seq_total)
names(cop_olig) <- c('abundances', 'seq_total')
cop_olig$abundances <- as.matrix(cop_olig$abundances)
cop_olig$rel.abundances <- cop_olig$abundances / cop_olig$seq_total
# saveRDS(cop_olig, NEON_cop_olig_abundances.path)


########## 3. assign nitrogen-cycling groups. ############

#### Set up N-cycle dataset from Albright 2008 ####

# Remove all columns except for taxonomy, environment, and pathways
N_cyclers <- N_cyclers_raw[,-c(1:2, 4, 11:14, 16:17)]
# Rename some pathways
setnames(N_cyclers,
         c("Nitrogen Fixation", "Assimilatory Nitrite to ammonia",
           "Dissimilatory Nitrite to Ammonia", "Assimilatory Nitrate to Nitrite",
           "Dissimilatory Nitrate to Nitrite"),
         c("N_fixation", "Assim_nitrite_reduction", "Dissim_nitrite_reduction",
           "Assim_nitrate_reduction", "Dissim_nitrate_reduction"))

# Treat "incomplete" pathways as if they are absent.
N_cyclers[N_cyclers == "complete"] <- 1
N_cyclers[N_cyclers == "incomplete" | N_cyclers == "None"] <- 0

# Grouping partial nitrification pathway with nitrification,
# and partial denitrification with denitrification.
N_cyclers[N_cyclers$Partial_Nitrification == 1,]$Nitrification <- 1
N_cyclers[N_cyclers$Partial_NO == 1 | N_cyclers$Partial_N2O == 1 |
            N_cyclers$Partial_N2 == 1,]$Denitrification <- 1

# Now we can remove the "partial" columns.
N_cyclers[, c("Partial_Nitrification",
              "Partial_NO",
              "Partial_N2O",
              "Partial_N2")] <- NULL

#### Assign taxa to functional groups ####

# create pathway columns
tax <- tax_save
pathway_names <- colnames(N_cyclers[9:15])
tax[, pathway_names] <- NA

# taxon assignments
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  
  # Classifications from literature search (multiple taxon levels)
  # if taxon is in classification data, assign it a present pathway
  has_pathway <- fg[fg$Classification == p,]$Taxon
  levels <- c("phylum", "class", "order", "family", "genus")
  for (j in 1:length(levels)) {
    k <- levels[j]
    if (nrow(tax[tax[[k]] %in% has_pathway,]) != 0) { 
      tax[tax[[k]] %in% has_pathway,][p] <- 1
    }
  }
  # genus + species must match any full species name
  if (nrow(tax[which(paste(tax$genus, tax$species) %in% has_pathway),]) != 0) {
    tax[which(paste(tax$genus, tax$species) %in% has_pathway),][p] <- 1
  }
  # Classifications from Albright et al. 2018 dataset (Genus-level only; reduced from species-level)
  has_pathway <- N_cyclers[N_cyclers[p] == 1,]$Genus # one species with pathway is enough to classify genus
  if (nrow(tax[which(tax$genus %in% has_pathway),][p]) != 0) {
    tax[which(tax$genus %in% has_pathway),][p] <- 1
  }
}

#Get seq abundances of each pathway
all_N_pathways <- list()
k <- data.table(cbind(tax, t(otu)))
for (i in 1:length(pathway_names)) {
  pathways <- list()
  n <- pathway_names[i]
  z <- k[k[[n]] == 1,]
  start <- ncol(tax) + 1
  out <- colSums(z[, start:ncol(z)])
  pathways[[i]] <- out
  pathways <- t(do.call('rbind', pathways))
  colnames(pathways) <- pathway_names[i]
  seq_total <- colSums(k[, start:ncol(k)])
  other <- seq_total - rowSums(pathways)
  pathways <- cbind(other, pathways)
  pathways <- list(pathways, seq_total)
  names(pathways) <- c('abundances', 'seq_total')
  pathways$rel.abundances <-
    pathways$abundances / pathways$seq_total
  all_N_pathways[[i]] <- pathways
}
#saveRDS(allpathways, NEON_N_cyclers_abundances.path)


########## 4. assign C-cycling groups. ###########

#### format dataset of cellulolytic taxa from Berlemont et al. ####

rownames(cell) <- cell$Strain
# get taxa with any enzymes for cellulolysis
cellulolytic <- cell[, c(5:7, 9:13)] #subset to names and pathways
cellulolytic$is.cell <- 0
cellulolytic[apply(cellulolytic, 1, function(x)
  any(x == 1)),]$is.cell <- 1
cellulolytic$genus <- word(rownames(cellulolytic), 1) # grab first word (genus)
cellulolytic[cellulolytic$genus == "Candidatus",]$genus <- # if first word is candidatus, grab two words
  word(rownames(cellulolytic[cellulolytic$genus == "Candidatus",]), 1, 2) 

# taxon assignments
levels <- c("phylum", "class", "order", "family", "genus")
for (i in 1:length(levels)) {
  p <- levels[i]
  tax[which(tax[[p]] %in% copiotrophs),]$group <- "copiotroph"
  tax[which(tax[[p]] %in% oligotrophs),]$group <- "oligotroph"
}


#### assign carbon-cycling taxa ####

tax <- tax_save
pathway_names <-
  c("Cellulolytic", "Chitinolytic", "Lignolytic", "Methanotroph") 
tax[, pathway_names] <- NA

# check if sample genus is in classification data, and that a classified pathway is present;
# assign those genera a present pathway
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  
  # Classifications from literature search (multiple taxon levels)
  has_pathway <- fg[which(fg$Classification == p),]$Taxon
  levels <- c("phylum", "class", "order", "family", "genus")
  for (j in 1:length(levels)) {
    k <- levels[j]
    if (nrow(tax[tax[[k]] %in% has_pathway,]) != 0) { 
      tax[tax[[k]] %in% has_pathway,][p] <- 1
    }
  }
  
  # genus + species must match any full species name
  if (nrow(tax[which(paste(tax$genus, tax$species) %in% has_pathway),]) != 0) {
    tax[which(paste(tax$genus, tax$species) %in% has_pathway),][p] <- 1
  }
  # Classifications from Berlemont et al. 2018 dataset (Genus-level only; reduced from species-level)
  if (p == "Cellulolytic") {
    has_pathway <- cellulolytic[cellulolytic$is.cell == 1,]$genus # one species with pathway is enough to classify genus
    if (nrow(tax[which(tax$genus %in% has_pathway),][p]) != 0) {
      tax[which(tax$genus %in% has_pathway),][p] <- 1
    }
  } # close cellulolytic section
}

# # check how many ended up without classifications.
# tax_classified <- tax[, 8:11]
# no_pathways <- tax_classified[apply(tax_classified, 1, function(x) !any(x == 1)),]
# nrow(no_pathways)
# [1] 151048 

#Get seq abundances of each pathway
all_C_pathways <- list()
k <- data.table(cbind(tax, t(otu)))
for (i in 1:length(pathway_names)) {
  pathways <- list()
  n <- pathway_names[i]
  z <- k[k[[n]] == 1,]
  start <- ncol(tax) + 1
  out <- colSums(z[, start:ncol(z)])
  pathways[[i]] <- out
  pathways <- t(do.call('rbind', pathways))
  colnames(pathways) <- pathway_names[i]
  seq_total <- colSums(k[, start:ncol(k)])
  other <- seq_total - rowSums(pathways)
  pathways <- cbind(other, pathways)
  pathways <- list(pathways, seq_total)
  names(pathways) <- c('abundances', 'seq_total')
  pathways$rel.abundances <-
    pathways$abundances / pathways$seq_total
  all_C_pathways[[i]] <- pathways
}
#saveRDS(allpathways, NEON_C_cyclers_abundances.path)

# combine each list
fg_abundances <- c(all_N_pathways, all_C_pathways, list(cop_olig))
fg_names <- list()
for (i in 1:12) {
  fg_names[[i]] <- colnames(fg_abundances[[i]][[1]])[2]
}
fg_names[[12]] <- "Cop_olig" #Cop_olig has one more column than the other 11 
names(fg_abundances) <- fg_names

# save all functional group abundances.
saveRDS(fg_abundances, NEON_fg_abundances_16S.path)
