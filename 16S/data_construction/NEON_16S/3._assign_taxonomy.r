#1. Assign taxonomy using dada2 in parallel. 2. Create table of most abundant genera.
#This script assumes you have a taxonomy table where:
#1. the row names are sample names.
#2. the column names are the actual unique sequences.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
library(doParallel)

#Here i load an OTU table with column names as unique sequences to assign.
otu <- readRDS(NEON_dada2_SV_table.path)
to_assign <- colnames(otu) #grab sequences to assign.

#specify output path here.
tax_output_path <- NEON_dada2_tax_table.path

#Everything from here below *should* just run and save where you told it to.
#download greengenes reference database.
cat('downloading green genes...\n')
greengenes.path <- paste0(data.dir,'gg_13_8_train_set_97.fa')
gg.download.link <- 'https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1'
cmd <- paste0('curl ',gg.download.link,' > ',greengenes.path)
system(cmd)
cat('greengenes download complete.\n')

#assign taxonomy.
tic()
cat('Assigning taxonomy using the RDP Classifier...\n')
out <- dada2::assignTaxonomy(to_assign,greengenes.path,multithread = T)
cat('Taxonomy assignment complete. ')
toc()

#save output as your taxonomy file.
saveRDS(out, tax_output_path)




tax <- readRDS(NEON_dada2_tax_table.path)
otu <- readRDS(NEON_dada2_SV_table.path)

# remove leading "k__" in taxonomy.
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

# for column names to be lower case.
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to bacteria or archaea from tax and otu table.
tax <- tax[!is.na(tax$kingdom),]
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

# normalize the otu table 
otu <- t(otu)
pro.function <- function(otu){
  for(i in 1:ncol(otu)){
    otu[,i] <- otu[,i] / sum(otu[,i])
  }
  return(otu)
}
otu <- pro.function(otu)

# make sure column sums are 1
colSums(otu)

# aggregate important classes and genera
# get most abundant genera in dataset.
genera <- unique(tax$genus)
test <- data.table(cbind(tax, otu))
seq.out <- list()
for(i in 1:length(genera)){
  z <- test[genus == genera[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  seq.out[[i]] <- out
}

# Count genus level abundance, grab some number of most abundant genera
seq.out <- do.call('rbind',seq.out)
counts <- rowSums(seq.out)
genera <- as.character(genera)
k <- data.table(cbind(genera,counts))
k$counts <- as.numeric(as.character(k$counts))
k <- k[order(-counts),]
k <- k[genera!=""&!is.na(genera),] #remove NA and empty genera
#grab genera of interest.
of_interest <- k$genera[1:n.gen]

# Get relative abundances of the most abundant genera
gen.list <- list()
for(i in 1:length(of_interest)){
  z <- data.table(cbind(tax,otu))
  z <- z[genus %in% of_interest[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  gen.list[[i]] <- out
}
gen.list <- data.frame(t(do.call('rbind',gen.list)))
colnames(gen.list) <- of_interest
gen.list$Mapping.ID <- rownames(gen.list)

saveRDS(gen.list, NEON_gen_abundances.path)
