# script to subset out only taxa that are found in 95%+ of samples
rm(list=ls())
source('paths.r')

#load the relative taxonomic data
data.list <- readRDS(emp_phylo.level.list_esv.path)

# find number for majority of the samples
#n_half <- ncol(data.list[[1]])/2 #same for every level of taxonomy
n_maj <- ncol(data.list[[1]]) * .95

#loop through each data set to subset out only taxa found in 50%+ of samples
for(i in 1:length(data.list)){
  count <- 0 
  num <- 1
  names <- character()
  #loop through and count number of samples taxa are in
  for(j in 1:nrow(data.list[[i]])){ #rows are taxa
    for(k in 1:ncol(data.list[[i]])){ #columns are samples
      if(data.list[[i]][j,k] > 0){
        count = count + 1
      }
    }
    if(count >= n_maj){
      names[num] <- rownames(data.list[[i]])[j]
      num = num + 1
    }
    count = 0
  }
  data.list[[i]] <- data.list[[i]][names,]
}

#save as RDS
saveRDS(data.list, emp_phylo.level.list_esv.path)
