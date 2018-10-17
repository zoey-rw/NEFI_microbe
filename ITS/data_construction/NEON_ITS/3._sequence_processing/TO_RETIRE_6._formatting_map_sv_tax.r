#hierarchical means core-plot-site.
#load SV table and map.
rm(list=ls())
library(runjags)
library(data.table)
source('paths.r')
source('NEFI_functions/fg_assign.r')
source('NEFI_functions/fg_assign_parallel.r')

#set output paths.----
output.path <- NEON_taxa_fg.path

#load data.----
sv <- readRDS(NEON_SV.table.path)
map <- readRDS(hierarch_filled.path)
map <- map$core.obs
fun <- readRDS(NEON_fun.path)

#subset to sv rows in map file. deprecatedVialID links map to sv.table.----
rownames(sv) <- substr(rownames(sv),4,nchar(rownames(sv)))
sv.id <- rownames(sv)
map <- map[map$deprecatedVialID %in% sv.id,]
 sv <-  sv[rownames(sv) %in% map$deprecatedVialID,]
rownames(map) <- map$deprecatedVialID

#put in same order
map <- map[order(rownames(map)),]
 sv <-  sv[order(rownames( sv)),]

#remove samples with less than 1000 reads from map and sv.----
map$seq.depth <- rowSums(sv)
map <- map[map$seq.depth > 1000,]
 sv <- sv[rownames(sv) %in% rownames(map),]

#kill SVs that no longer have any sequences or are now singletons.----
sv.filter <- data.frame(colnames(sv),colSums(sv))
colnames(sv.filter) <- c('sv.ID','n.seqs')
to_remove <- sv.filter[sv.filter$n.seqs < 2,]$sv.ID
 sv <- sv[,!(colnames(sv) %in% to_remove)]
fun <- fun[!(rownames(fun) %in% to_remove),]

#remove taxa that do not assign to fungi.----
to_remove <- rownames(fun[is.na(fun$kingdom),])
fun <- fun[!(rownames(fun) %in% to_remove),]
 sv <-  sv[,!(colnames(sv) %in% to_remove) ]
 
#Assign functional groups.----
#Colin double checked everything stays in order.
#Things can't be both ECM and SAP in dirlichet. Creates a sum to 1 problem.
tax <- data.table(fun)
tax[grep('Ectomycorrhizal', tax$guild),Ectomycorrhizal := 1]
tax[is.na(Ectomycorrhizal),            Ectomycorrhizal := 0]
tax[grep('Saprotroph', tax$guild),          Saprotroph := 1]
tax[is.na(Saprotroph),                      Saprotroph := 0]
tax[grep('Arbuscular', tax$guild),          Arbuscular := 1]
tax[is.na(Arbuscular),                      Arbuscular := 0]
tax[grep('Patho', tax$guild),                 Pathogen := 1]
tax[is.na(Pathogen),                          Pathogen := 0]
#If you are ECTO you can't be SAP
tax[Ectomycorrhizal == 1, Saprotroph := 0]
#if you are ecto or sap you not path.
tax[Ectomycorrhizal == 1,   Pathogen := 0]
tax[     Saprotroph == 1,   Pathogen := 0]
#convert back to data.frame.
tax <- as.data.frame(tax)
rownames(tax) <- rownames(fun)

#Get counts of fucntional groups. (ECM, AM, SAP, PATH)
function_groups <- c('Ectomycorrhizal','Arbuscular','Saprotroph','Pathogen')
fun.list <- list()
tax$map.ID <- rownames(tax)
z <- data.table(cbind(tax,t(sv)))
for(i in 1:length(function_groups)){
  k <- z[eval(parse(text=function_groups[i]))== 1,]
  start <- ncol(tax) + 1
  out <- colSums(k[,start:ncol(k)])
  fun.list[[i]] <- out
}
fun.list <- data.frame(t(do.call('rbind',fun.list)))
colnames(fun.list) <- function_groups


#get Cosmopolitan genera presence/absence and counts.----
genera <- unique(tax$genus)
test <- data.table(cbind(tax, t(sv)))
abundance <- list()
 presence <- list()
for(i in 1:length(genera)){
  z <- test[genus == genera[i],]
  start <- ncol(tax) + 1
  out1 <- colSums(z[,start:ncol(z)])
  out2 <- ifelse(out1 > 0, 1, 0)
  abundance[[i]] <- out1
   presence[[i]] <- out2
}
#Count genus level abundance, grab some number of most abundant genera
count.out <- do.call('rbind',abundance)
 pres.out <- do.call('rbind',presence)
rownames(count.out) <- genera
rownames( pres.out) <- genera
colnames(count.out) <- rownames(sv)
colnames( pres.out) <- rownames(sv)

#get most cosmopolitan genera in at least half of samples.
N.samps <- nrow(sv)
cosmo <- rowSums(pres.out)
names(cosmo) <- NULL
cosmo <- data.frame(cosmo)
cosmo <- data.frame(genera,cosmo)
cosmo <- cosmo[order(cosmo$cosmo, decreasing = T),]
#subset to genera present in at least half of samples.
cosmo <- cosmo[cosmo$cosmo > N.samps/2,]
cosmo.counts <- count.out[rownames(count.out) %in% cosmo$genera,]
cosmo.counts <- data.frame(t(cosmo.counts))


#merge in cosmopolitan genera counts.----
fun.list$ID <- rownames(fun.list)
cosmo.counts$ID <- rownames(cosmo.counts)
counts <- merge(fun.list, cosmo.counts)
counts$seq.depth <- rowSums(sv)
rownames(counts) <- counts$ID

#calculate relative counts.----
rel.counts <- counts
rel.counts[,2:(ncol(rel.counts) - 1)] <- rel.counts[,2:(ncol(rel.counts) - 1)] / rel.counts$seq.depth
rel.counts[is.na(rel.counts)] <- 0

#get other IDs in here. deprecatedVialID does not match all products.----
    counts <- merge(counts,map[,c('deprecatedVialID','geneticSampleID')], by.x = 'ID', by.y = 'deprecatedVialID', all.x = T)
rel.counts <- merge(rel.counts, map[,c('deprecatedVialID','geneticSampleID')], by.x = 'ID', by.y = 'deprecatedVialID', all.x=T)
  
#save output.----
dat.out <- list(counts,rel.counts)
names(dat.out) <- c('counts','rel.counts')
saveRDS(dat.out, output.path)
