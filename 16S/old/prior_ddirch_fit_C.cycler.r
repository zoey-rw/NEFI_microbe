# Fit dirlichet models to C-cycling groups of bacteria/archaea from Bahram et al. 
# Model fits separately for each of 4 C-cycling groups, since they are not mutually-exclusive.
# Taxa were categorized through a lit review.
# No hierarchy required, as everything is observed at the site level. Each observation is a unique site.

#clear environment
rm(list = ls())
library(data.table)
library(runjags)
library(doParallel)
library(future)
library(furrr) 
library(magrittr)
plan(multiprocess)
source('paths.r')
source('NEFI_functions/ddirch_site.level_JAGS.r')
source('NEFI_functions/ddirch_site.level_JAGS_int.only.r')
source('NEFI_functions/crib_fun.r')


#detect and register cores.
n.cores <- detectCores()
registerDoParallel(cores=n.cores)

#set output path.----
output.path <- bahram_16S_prior_C_cycle_JAGSfits 

# load Bahram metadata.----
m <- data.table(readRDS(bahram_metadata.path))
setnames(m, "PH", "pH")

# load Bahram abundance data
a <- readRDS(prior_C_cyclers_abundances.path)

# load covariate selection data.
covs <- readRDS(bahram_16S_prior_fg_cov.selection_JAGS)

# set pathway names
pathway_names <- list()
for (i in 1:4) {
  pathway_names[[i]] <- colnames(a[[i]]$abundances)[2]
}

all_C_cycler_models <- vector("list", 4)

# loop through each C_cycle group.
for (i in 1:length(a)) {
  y <- a[[i]]
  y <- y$abundances
  
  d <- m[,.(Run,pC,cn,pH,Ca,Mg,P,K,pN,moisture,NPP,map,mat,forest,conifer,relEM)] 
  d <- d[complete.cases(d),] #optional. This works with missing data.
  y <- y[rownames(y) %in% d$Run,]
  
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
  x <- d
  rownames(x) <- x$Run
  x$Run <- NULL
  intercept <- rep(1, nrow(x))
  x <- cbind(intercept, x)
  
  #IMPORTANT: LOG TRANSFORM MAP.
  #log transform map, magnitudes in 100s-1000s break JAGS code.
  x$map <- log(x$map)
  
  # define multiple subsets
  # get covariates from model selection output
  covariates <- c("intercept",rownames(covs[[2]][[i]])) # C_cyclers are second item; each of 4 pathways has different covariates
  cols <- which(colnames(x) %in% covariates)
  x.cov_select <- x[,cols, with=FALSE]
  x.all  <- x[,.(intercept,pC,cn,pH,Ca,Mg,P,K,pN,moisture,NPP,mat,map,forest,conifer,relEM)] # all nutrients, no moisture
  x.list <- list(x.cov_select,x.all)
  
  #fit model using function.
  #This take a long time to run, probably because there is so much going on.
  #fit <- site.level_dirlichet_jags(y=y,x_mu=x,adapt = 50, burnin = 50, sample = 100)
  #for running production fit on remote.
  output.list<-
    (1:length(x.list)) %>%	
    future_map(function(i){
      fit <- site.level_dirlichet_jags(y=y,x_mu=x.list[i],adapt = 200, burnin = 1000, sample = 1000, parallel = T)
      return(fit)
    })
  
  #get intercept only fit.
  output.list[[length(x.list) + 1]] <- site.level_dirlichet_intercept.only_jags(y=y, silent.jags = T)
  
  #name the items in the list
  names(output.list) <- c('cov_select','all.preds','int.only')
  all_C_cycler_models[[i]] <- output.list
  cat(paste0('C-cycler fit completed for ', pathway_names[[i]]))
}

cat('Saving fit...\n')
saveRDS(all_C_cycler_models, output.path)
cat('Script complete. \n')