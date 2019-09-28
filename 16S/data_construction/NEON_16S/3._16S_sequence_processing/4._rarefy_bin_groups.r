#NEON taxonomic group tables for all phylogenetic levels - 16S.
#common = found in greater than 50% of prior samples.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')
source('NEFI_functions/common_group_quantification.r')

#set output path.----
output.path <- NEON_16S_phylo_fg_abundances.path
# load NEON SV table as otu file
otu <- readRDS(NEON_dada2_SV_table.path)
# load NEON taxonomy
tax <- readRDS(NEON_dada2_tax_table.path)
# load reference taxonomy-to-function file.
tax_fun <- readRDS(paste0(pecan_gen_16S_dir, "reference_data/bacteria_tax_to_function.rds"))
# load reference common phyla from prior.----
ref <- readRDS(delgado_16S_common_phylo_fg_abun.path)

########## 1. prep taxonomic data. ############

# remove leading "k__" in taxonomy.
for (i in 1:ncol(tax)) {
  tax[, i] <- substring(tax[, i], 4)
}

# for column names to be lower case.
tax <- as.data.frame(tax)
#colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$Kingdom == 'Bacteria',] # removes ~900 counts
otu <- otu[, colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

# rarefy otu table
set.seed(5) # so that rarefaction is repeatable.
otu <- otu[rowSums(otu) >= 10000,]
otu <- vegan::rrarefy(otu, 10000)

# assign function to taxonomy
pathway_names <- colnames(tax_fun)[3:15]
tax[, pathway_names] <- "other"
# taxon assignments
for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  
  # Classifications from literature search (multiple taxon levels)
  # I'm so sorry for anyone observing this nested for-loop in the future
  has_pathway <- tax_fun[tax_fun[,p] == 1,]
  levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (j in 1:length(levels)) {
    taxon_level <- levels[j]
    has_pathway_taxon_level <- has_pathway[has_pathway$Taxonomic.level==taxon_level,]
    if (taxon_level == "Species") {
      if(nrow(tax[which(paste(tax$Genus, tax$Species) %in% has_pathway_taxon_level$Taxon),]) > 0) {
        tax[which(paste(tax$Genus, tax$Species) %in% has_pathway_taxon_level$Taxon),][,p] <- p
      }
    } else {
      if (nrow(tax[tax[[taxon_level]] %in% has_pathway_taxon_level$Taxon,]) > 0){
        tax[tax[[taxon_level]] %in% has_pathway_taxon_level$Taxon,][,p] <- p
      }
    }
  }
}



#check that things are properly ordered.----
if(sum(rownames(tax) == colnames(otu)) != ncol(otu)){
  cat('Stop. rownames of tax-fun table do not match column names of SV table.\n')
}

#get each level of taxonomy output.----
of_interest <- colnames(tax)[!colnames(tax) %in% c("Kingdom","Species")]
all_taxa_out <- list()
for(i in 1:length(of_interest)){
  # reference <- data.table(ref[[i]]$group_frequencies)
  # reference <- as.character(reference[sample_frequency > 0.95]$groups)
  reference <- data.table(ref[[i]]$abundances)
  reference <- colnames(reference)[!colnames(reference) %in% "other"]
  sv = otu
tax = tax
groups = reference
tax_level = of_interest[i]
ref_filter = T
all_taxa_out[[i]] <- common_group_quantification(otu,
                                                   tax,
                                                   reference,
                                                   of_interest[i],
                                                   ref_filter = T)
}
names(all_taxa_out) <- of_interest

#save output.----
saveRDS(all_taxa_out,output.path) 
