#aggregate tedersoo functional groups.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/fg_assign.r')
library(data.table)

#set output paths.----
   fg.output.path <- tedersoo_ITS_fg_list.path
hydro.output.path <- tedersoo_ITS_hydro_list.path
yeast.output.path <- tedersoo_ITS_yeast_list.path

#load data.----
map <- readRDS(tedersoo_ITS_clean_map.path)
otu <- readRDS(ted_2014_SV.table.path)
tax <- readRDS(ted_2014_tax.path)
em.trait <- data.table(read.csv(em_traits.path)) #hydrophillic/phobic assignments from Hobbie.

#subset otu file to match mapping file.----
otu <- otu[rownames(otu) %in% map$SRR.id,]

#format taxonomy table.----
#remove leading "k__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

#for column names to be lower case.
colnames(tax) <- tolower(colnames(tax))

#remove taxa that do not assign to fungi from tax and otu table.
tax <- tax[tax$kingdom == 'Fungi',]
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

#Assign functional groups using FUNGuild, append to taxonomy table.----
tax <- fg_assign(tax)

#Categorize tax table as ECM, SAP, Arbuscular, saprotroph.----
tax <- data.table(tax) #data.table makes this go fast.
tax[grep('Ectomycorrhizal', tax$guild),Ectomycorrhizal := 1]
tax[is.na(Ectomycorrhizal),            Ectomycorrhizal := 0]
tax[grep('Saprotroph', tax$guild),          Saprotroph := 1]
tax[is.na(Saprotroph),                      Saprotroph := 0]
tax[grep('Arbuscular', tax$guild),          Arbuscular := 1]
tax[is.na(Arbuscular),                      Arbuscular := 0]
tax[grep('Patho', tax$guild),                 Pathogen := 1]
tax[is.na(Pathogen),                          Pathogen := 0]
tax[grep('Yeast', tax$growthForm),               Yeast := 1]
tax[grep('Facultative Yeast', tax$growthForm),   Yeast := 0]
tax[grep('Facultative Yeast', tax$growthForm),   Facultative_Yeast := 1]
tax[is.na(Facultative_Yeast),         Facultative_Yeast:= 0]
tax[is.na(Yeast),                                Yeast := 0]

#Things can't have multiple functional assignments in dirlichet. Creates a sum to 1 problem.
#If you are ECTO you can't be SAP
tax[Ectomycorrhizal == 1, Saprotroph := 0]
#if you are ecto or sap you not path.
tax[Ectomycorrhizal == 1,   Pathogen := 0]
tax[     Saprotroph == 1,   Pathogen := 0]

#categorize tax table hydrophobic/hydrophillic status.----
tax <- data.table(tax)
tax[genus %in% em.trait[ hydrophobic == 1,]$genus, hydrophobic := 1]
tax[genus %in% em.trait[hydrophillic == 1,]$genus,hydrophillic := 1]
tax[is.na( hydrophobic),  hydrophobic := 0]
tax[is.na(hydrophillic), hydrophillic := 0]

#Get abundances and relative abundances of hydrophillic/phobic groups.----
#hydrophillic/phobic.
z <- data.table(cbind(tax,t(otu)))
phil <- z[hydrophillic == 1,]
phob <- z[hydrophobic  == 1,]
start <- ncol(tax) + 1
out.phil <- colSums(phil[,start:ncol(phil)])
out.phob <- colSums(phob[,start:ncol(phob)])
hydro.out <- data.frame(t(rbind(out.phil,out.phob)))
colnames(hydro.out) <- c('hydrophillic','hydrophobic')
seq_total <- colSums(z[,start:ncol(z)])
other <- seq_total - rowSums(hydro.out)
hydro.out <- cbind(other,hydro.out)
hydro.list <- list(hydro.out, seq_total)
names(hydro.list) <- c('abundances','seq_total')
hydro.list$rel.abundances <- hydro.list$abundances / hydro.list$seq_total

#get abundances and relative abundances of functional groups (ECM, AM, SAP, WR).----
function_groups <- c('Ectomycorrhizal','Arbuscular','Saprotroph','Pathogen')
fun.list <- list()
k <- data.table(cbind(tax,t(otu)))
for(i in 1:length(function_groups)){
  z <- k[eval(parse(text=function_groups[i]))== 1,]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  fun.list[[i]] <- out
}
fun.list <- data.frame(t(do.call('rbind',fun.list)))
colnames(fun.list) <- function_groups
seq_total <- colSums(k[,start:ncol(k)])
other <- seq_total - rowSums(fun.list)
fun.list <- cbind(other,fun.list)
fun.list <- list(fun.list,seq_total)
names(fun.list) <- c('abundances','seq_total')
fun.list$rel.abundances <- fun.list$abundances / fun.list$seq_total
#fun.list$Mapping.ID <- rownames(fun.list)

#get abundances of yeast and non-yeast.----
function_groups <- c('Yeast','Facultative_Yeast')
yeast.list <- list()
for(i in 1:length(function_groups)){
  z <- k[eval(parse(text=function_groups[i]))== 1,]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  yeast.list[[i]] <- out
}
yeast.list <- do.call(cbind,yeast.list)
colnames(yeast.list) <- function_groups
yeast.abundances <- cbind(seq_total - rowSums(yeast.list),yeast.list)
colnames(yeast.abundances) <- c('other',function_groups)
yeast.rel.abundances <- yeast.abundances / seq_total
yeast.list <- list(yeast.abundances, yeast.rel.abundances,seq_total)
names(yeast.list) <- c('abundances','rel.abundances','seq_total')

#save outputs.----
saveRDS(hydro.list, hydro.output.path)
saveRDS(  fun.list,    fg.output.path)
saveRDS(yeast.list, yeast.output.path)
cat('YOU DID IT. GREAT JOB. <3 Colin.')

