#hierarchical means core-plot-site.
#load SV table and map.
rm(list=ls())
library(runjags)
source('paths.r')
source('NEFI_functions/fg_assign.r')
source('NEFI_functions/fg_assign_parallel.r')
sv <- readRDS(NEON_SV.table.path)
map <- readRDS(hierarch_filled.path)
map <- map$core.obs
tax <- readRDS(NEON_tax.path)

#subset to sv rows in map file. deprecatedVialID links map to sv.table
rownames(sv) <- substr(rownames(sv),4,nchar(rownames(sv)))
sv.id <- rownames(sv)
map <- map[map$deprecatedVialID %in% sv.id,]
 sv <-  sv[rownames(sv) %in% map$deprecatedVialID,]

#remove taxa that do not assign to fungi.
to_remove <- rownames(tax[is.na(tax$Kingdom),])
tax <- tax[!(rownames(tax) %in% to_remove),]
 sv <-  sv[,!(colnames(sv) %in% to_remove) ] 
 
#assign functional groups.
colnames(tax) <- tolower(colnames(tax))
#remove leading "k__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}
tax.test <- tax[1:1000,]

#TESTING ASSIGNMENT
tax <- readRDS(NEON_tax.path)
to_remove <- rownames(tax[is.na(tax$Kingdom),])
tax <- tax[!(rownames(tax) %in% to_remove),]
test <- fg_assign_parallel(tax.test, n.cores = 2)

