# clear environment and set working directory
rm(list=ls())
source('paths.r')
source('NEFI_functions/pro_norm_otu.r')

#load the taxonomic assignments
tax <- readRDS(emp_tax.path)
esv <- readRDS(emp_esv_clean.path)
map <- readRDS(emp_map_clean.path)

#remove the unwanted string characters
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4) # removes first three characters (e.g. "p__" in front of phylum assignment)
}
for(i in 1:ncol(tax)){
  tax[,i] <- gsub("\\[|\\]","",tax[,i]) # removes square brackets around assignments
}

#make taxonomy vectors
phylum <- tax$Phylum[match(rownames(esv),rownames(tax))]
 class <-  tax$Class[match(rownames(esv),rownames(tax))]
 order <-  tax$Order[match(rownames(esv),rownames(tax))]
family <- tax$Family[match(rownames(esv),rownames(tax))]
 genus <-  tax$Genus[match(rownames(esv),rownames(tax))]

#aggregate esv file by taxonomy
phylum_esv <- aggregate(esv,by=list(phylum),FUN=sum)
 class_esv <- aggregate(esv,by=list(class) ,FUN=sum)
 order_esv <- aggregate(esv,by=list(order) ,FUN=sum)
family_esv <- aggregate(esv,by=list(family),FUN=sum)
 genus_esv <- aggregate(esv,by=list(genus) ,FUN=sum)

#make taxonomy the rownames of new esv file
rownames(phylum_esv) <- phylum_esv$Group.1; phylum_esv$Group.1 <- NULL
rownames( class_esv) <-  class_esv$Group.1;  class_esv$Group.1 <- NULL
rownames( order_esv) <-  order_esv$Group.1;  order_esv$Group.1 <- NULL
rownames(family_esv) <- family_esv$Group.1; family_esv$Group.1 <- NULL
rownames( genus_esv) <-  genus_esv$Group.1;  genus_esv$Group.1 <- NULL

# select only samples wanted
phylum_esv <- phylum_esv[colnames(phylum_esv) %in% map$`#SampleID`]
 class_esv <-  class_esv[colnames(class_esv ) %in% map$`#SampleID`]
 order_esv <-  order_esv[colnames(order_esv ) %in% map$`#SampleID`]
family_esv <- family_esv[colnames(family_esv) %in% map$`#SampleID`]
 genus_esv <-  genus_esv[colnames(genus_esv ) %in% map$`#SampleID`]
 

data.list <- list(phylum_esv,class_esv,order_esv,family_esv,genus_esv)
names(data.list) <- c('phylum','class','order','family','genus')
#proportionally normalize each table.
for(i in 1:length(data.list)){
  data.list[[i]] <- pro_norm_otu(data.list[[i]])
}

#save esv files
saveRDS(data.list, emp_phylo.level.list_esv.path) 
