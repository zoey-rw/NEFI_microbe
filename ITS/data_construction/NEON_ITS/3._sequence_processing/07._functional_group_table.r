#NEON functional group tables.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')

#set output paths.----
output.path <- NEON_taxa_fg.path

#load data.----
sv <- readRDS(NEON_SV.table_clean.path)
map <- readRDS(hierarch_filled.path)
map <- map$core.obs
fun <- readRDS(NEON_fun_clean.path)

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

#Get together abundances and seq_total.----
seq_total <- colSums(z[,start:ncol(z)])
seq_total <- seq_total[order(match(names(seq_total),rownames(fun.list)))]
other <- seq_total - rowSums(fun.list)
abundances <- cbind(other,fun.list)
rel.abundances <- abundances / seq_total

#get other IDs in here. deprecatedVialID does not match all products.----
    abundances$deprecatedVialID <- rownames(    abundances)
rel.abundances$deprecatedVialID <- rownames(rel.abundances)

    abundances <- merge(    abundances, map[,c('deprecatedVialID','geneticSampleID')], 
                            by.x = 'deprecatedVialID', by.y = 'deprecatedVialID', all.x = T)
rel.abundances <- merge(rel.abundances, map[,c('deprecatedVialID','geneticSampleID')], 
                            by.x = 'deprecatedVialID', by.y = 'deprecatedVialID', all.x = T)

#save output.----
dat.out <- list(abundances,rel.abundances,seq_total)
names(dat.out) <- c('abundances','rel.abundances','seq_total')
saveRDS(dat.out, output.path)
