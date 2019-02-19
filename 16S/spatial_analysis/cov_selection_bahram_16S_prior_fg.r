#covariate selection on bahram 16Sprior functional groups.
#linear combination of as many predictors as we can think of, no interactions.
#clear environment, load paths and functions.
rm(list=ls())
source('NEFI_functions/covariate_selection_JAGS.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/tic_toc.r')
source('paths.r')

#load data.
d <- data.table::data.table(readRDS(bahram_metadata.path))
d <- d[,.(Run,pC,cn,pH,Ca,Mg,P,K,moisture,NPP,map,mat,forest,conifer,relEM)]
d <- d[complete.cases(d),] #optional. This works with missing data.

output.path <- bahram_16S_prior_fg_cov.selection_JAGS


# load Bahram functional group abundances
a1 <- readRDS(prior_N_cyclers_abundances.path)
a2 <- readRDS(prior_C_cyclers_abundances.path)
a3 <- list(readRDS(prior_cop_olig_abundances.path))

# combine these three lists of lists - just get abundances
a <- do.call(c, list(a1, a2, a3))
fg_all <- sapply(a, "[[", 1)

#fg_all <- c("N_cyclers", "C_cyclers", "Cop_olig")
all_fg_output <- list()

for (f in 1:length(fg_all)) {
  fg_output <- list()
  fg <- fg_all[f]
  
  #organize y data
  y <- fg
  y <- y[[1]]
  y <- y[rownames(y) %in% d$Run,]
  d <- d[d$Run %in% rownames(y)]
  #order abundance table to match the metadata file
  y <- y[match(d$Run, rownames(y)),]
  if(!sum(rownames(y) == d$Run) == nrow(y)){
    cat('Warning. x and y covariates not in the same order!')
  }
  
  #Get relative counts by adding 1 to all observations (can't handle zeros).----
  y <- y + 1
  y <- y/rowSums(y)
  y <- as.data.frame(y)
  
  #Drop in intercept, setup predictor matrix.
  d$intercept <- rep(1,nrow(d))
  d$map <- log(d$map)
  x <- as.data.frame(d[,.(intercept,pC,cn,pH,Ca,Mg,P,K,moisture,NPP,map,mat,forest,conifer,relEM)])
  x.no.nutr <- as.data.frame(d[,.(intercept,pC,cn,pH,moisture,NPP,map,mat,forest,conifer,relEM)])
  
  covs <- covariate_selection_JAGS(y=y,x_mu=x, n.adapt = 300, n.burnin = 1000, 
                                   n.sample = 1000,parallel = F)
  cat(paste0("First covariate selection complete for group: ", colnames(y)[2]))
  covs.no.nutr <- covariate_selection_JAGS(y=y,x_mu=x.no.nutr, n.adapt = 300, n.burnin = 1000, 
                                           n.sample = 1000, parallel = F)
  cat(paste0("Second (no.nutr) covariate selection complete for group: ", colnames(y)[2]))
  
  fg_output <- list(covs, covs.no.nutr)
  names(fg_output) <- c(paste0("Cov_select ", colnames(y)[2]), paste0("No.nutr cov_select ", colnames(y)[2]))
  cat(paste0("All covariate selection complete for group: ", colnames(y)[2]))
  all_fg_output[[f]] <- fg_output
  #return(fg_output)
  
} # end functional group loop


# set pathway names
group_names <- list()
for (i in 1:12) {
  group_names[[i]] <- colnames(fg_all[[i]])[2]
}
group_names[[12]] <- "Cop_olig" #Cop_olig has one more column than the other 11 
# add names to each item
names(all_fg_output) <- group_names
# names(covs) <- c("N_cycling","C_cycling","Cop_olig")
# names(covs[[1]]) <- c("Assim_nitrite_reduction", "Dissim_nitrite_reduction", "Assim_nitrate_reduction", 
#                       "N_fixation", "Dissim_nitrate_reduction", "Nitrification", "Denitrification")
# names(covs[[2]]) <- c("Cellulolytic", "Chitinolytic", "Lignolytic", "Methanotroph")
# names(covs[[3]]) <- "Cop_olig"

# save the output
saveRDS(all_fg_output, bahram_16S_prior_fg_cov.selection_JAGS)
