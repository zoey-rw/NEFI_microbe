#NEON functional group tables.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')

#set output paths.----
functional.output.path <- NEON_ITS_fastq_taxa_fg.path
     yeast.output.path <- NEON_ITS_fastq_yeast_taxa.path

#load data.----
sv <- readRDS(NEON_ITS_fastq_SV.table_clean.path)
map <- readRDS(hierarch_filled.path)
map <- map$core.obs
fun <- readRDS(NEON_ITS_fastq_fun_clean.path)

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
#assign yeast / non-yeast.
tax[grep('Yeast', tax$growthForm),               yeast := 1]
tax[is.na(yeast),                                yeast := 0]
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

#get a yeast list.
yeast.list <- list()
function_groups <- c('yeast')
z <- data.table(cbind(tax, t(sv)))
for(i in 1:length(function_groups)){
  k <- z[eval(parse(text=function_groups[i]))== 1,]
  start <- ncol(tax) + 1
  out <- colSums(k[,start:ncol(k)])
  yeast.list[[i]] <- out
}
yeast.list <- data.frame(t(do.call('rbind',yeast.list)))
colnames(yeast.list) <- function_groups


#Get together abundances and seq_total.----
seq_total <- colSums(z[,start:ncol(z)])
seq_total <- seq_total[order(match(names(seq_total),rownames(fun.list)))]
seq_total.yeast <- seq_total[order(match(names(seq_total), rownames(yeast.list)))]
other <- seq_total - rowSums(fun.list)
other.yeast <- seq_total.yeast - rowSums(yeast.list)
abundances <- cbind(other,fun.list)
yeast.abundances <- cbind(other.yeast, yeast.list)
rel.abundances <- abundances / seq_total

#get other IDs in here. deprecatedVialID does not match all products.----
     abundances$geneticSampleID <- rownames(    abundances)
 rel.abundances$geneticSampleID <- rownames(rel.abundances)
yeast.abundances$geneticSampleID <- rownames(yeast.abundances)

#save output.----
dat.out <- list(abundances,rel.abundances,seq_total)
names(dat.out) <- c('abundances','rel.abundances','seq_total')
yeast.out <- list(yeast.abundances, seq_total.yeast)
names(yeast.out) <- c('abundances','seq_total')
saveRDS(dat.out, functional.output.path)
saveRDS(yeast.out, yeast.output.path)
