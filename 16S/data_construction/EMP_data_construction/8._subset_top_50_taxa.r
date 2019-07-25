rm(list=ls())
source('paths.r')

#load the relative taxonomic data
data.list <- readRDS(emp_phylo.level.list_esv.path)
data.comp <- readRDS(emp_phylo.level.list_esv_comp.path)

for(i in 1:length(data.list)){
  data.comp[[i]] <- subset(data.comp[[i]], rownames(data.comp[[i]]) %in% rownames(data.list[[i]]))
  for(j in 1:nrow(data.comp[[i]])){
    data.comp[[i]]$Sum[j] <- sum(data.comp[[i]][j,])
  }
  #sort by total abundance
  most_abund <- data.comp[[i]][order(-data.comp[[i]]$Sum),]
  if(nrow(most_abund)>50){
    most_abund <- most_abund[1:50,]
  }
  #subset 50% data with the 50 most abundant
  data.list[[i]] <- data.list[[i]][rownames(most_abund),] 
}

#save objects as RDS
saveRDS(data.list, emp_phylo.level.list_esv.path)
